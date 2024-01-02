"""
Functions for fitting to data.

!!! warning
    This submodule is in beta. It's API may change.
"""
module Fit
using ..Korg, LineSearches, Optim
using Interpolations: LinearInterpolation
using ForwardDiff
using Trapz

# used by scale and unscale for some parameters
function tan_scale(p, lower, upper) 
    if !(lower <= p <= upper)
        throw(ArgumentError("p=$p is not in the range $lower to $upper"))
    end
    tan(π * (((p-lower)/(upper-lower))-0.5))
end
tan_unscale(p, lower, upper) = (atan(p)/π + 0.5)*(upper - lower) + lower

# these are the parmeters which are scaled by tan_scale
const tan_scale_params = Dict(
    "epsilon" => (0, 1),
    # we can't get these directly from Korg.get_atmosphere_archive() because it will fail in the 
    # test environent, but they are simply the boundaries of the SDSS marcs grid used by
    # Korg.interpolate_marcs.
    "Teff" => (2800, 8000),
    "logg" => (-0.5, 5.5),
    "m_H" => (-2.5, 1),
    # this allows all the atmospheres supported by the grid, but also many that are not.
    # alpha will be clamped to the nearest supported value.
    "alpha_H" => (-3.5, 2), 
    map(Korg.atomic_symbols) do el
        el => (-10, +4)
    end...
)

"""
Rescale each parameter so that it lives on (-∞, ∞).
"""
function scale(params::Dict)
    map(collect(params)) do (name, p)
        name => if name in keys(tan_scale_params)
            tan_scale(p, tan_scale_params[name]...)
        elseif name in ["vmic", "vsini"]
            tan_scale(sqrt(p), 0, sqrt(250))
        else
            @error "$name is not a parameter I know how to scale."
        end
    end |> Dict
end

"""
Unscale each parameter so that it lives on the appropriate range instead of (-∞, ∞).
"""
function unscale(params::Dict)
    map(collect(params)) do (name, p)
        name => if name in keys(tan_scale_params)
            tan_unscale(p, tan_scale_params[name]...)
        elseif name in ["vmic", "vsini"]
            tan_unscale(p, 0, sqrt(250))^2
        else
            @error "$name is not a parameter I know how to unscale."
        end
    end |> Dict
end

"""
Synthesize a spectrum, returning the flux, with LSF applied, resampled, and rectified.  This is 
used by fitting routines. See [`Korg.synthesize`](@ref) to synthesize spectra as a Korg user.
"""
function synthetic_spectrum(synthesis_wls, linelist, LSF_matrix, params;
                     line_buffer=10)
    specified_abundances = Dict([p for p in pairs(params) if p.first in Korg.atomic_symbols])
    alpha_H = "alpha_H" in keys(params) ? params["alpha_H"] : params["m_H"]
    A_X::Vector{valtype(params)} = Korg.format_A_X(params["m_H"], alpha_H, specified_abundances; solar_relative=true)

    # clamp_abundances clamps M_H, alpha_M, and C_M to be within the atm grid
    atm = Korg.interpolate_marcs(params["Teff"], params["logg"], A_X; clamp_abundances=true, perturb_at_grid_values=true)

    sol = Korg.synthesize(atm, linelist, A_X, synthesis_wls; vmic=params["vmic"], line_buffer=line_buffer, 
                          electron_number_density_warn_threshold=Inf)
    F = sol.flux ./ sol.cntm
    F = Korg.apply_rotation(F, synthesis_wls, params["vsini"], params["epsilon"])
    LSF_matrix * F
end

"""
Validate fitting parameters, and insert default values when needed. Used by [`fit_spectrum`](@ref).

these can be specified in either initial_guesses or fixed_params, but if they are not, these values
 are inserted into fixed_params
"""
function validate_params(initial_guesses::Dict, fixed_params::Dict;
                         required_params = ["Teff", "logg"],
                         default_params = Dict("m_H"=>0.0, "vsini"=>0.0, "vmic"=>1.0, "epsilon"=>0.6),
                         allowed_params = Set(["alpha_H" ; required_params ; keys(default_params)... ; Korg.atomic_symbols]))
    # convert all parameter values to Float64
    initial_guesses = Dict(string(p[1]) => Float64(p[2]) for p in pairs(initial_guesses))
    fixed_params = Dict(string(p[1]) => Float64(p[2]) for p in pairs(fixed_params))

    # check that all required params are specified
    all_params = keys(initial_guesses) ∪ keys(fixed_params)
    for param in required_params
        if !(param in all_params)
            throw(ArgumentError("Must specify $param in either starting_params or fixed_params. (Did you get the capitalization right?)"))
        end
    end

    # check that all the params are recognized
    unknown_params = filter!(all_params) do param
        param ∉ allowed_params
    end
    if length(unknown_params) > 0
        throw(ArgumentError("These parameters are not recognized: $(unknown_params)"))
    end 

    # filter out those that the user specified, and fill in the rest
    default_params = filter(collect(pairs(default_params))) do (k, v)
        !(k in keys(initial_guesses)) && !(k in keys(fixed_params))
    end |> Dict
    fixed_params = merge(default_params, fixed_params) 

    # check that no params are both fixed and initial guesses
    let keys_in_both = collect(keys(initial_guesses) ∩ keys(fixed_params))
        if length(keys_in_both) > 0
            throw(ArgumentError("These parameters: $(keys_in_both) are specified as both initial guesses and fixed params."))
        end 
    end

    initial_guesses, fixed_params
end

# make it possible to use dicts instead of NamedTuples for the python people
validate_params(initial_guesses::Dict, fixed_params::NamedTuple; kwargs...) = 
    validate_params(initial_guesses, _namedtuple_to_dict(fixed_params); kwargs...)
validate_params(initial_guesses::NamedTuple, fixed_params=Dict{String, Float64}(); kwargs...) = 
    validate_params(_namedtuple_to_dict(initial_guesses), fixed_params; kwargs...)

function _namedtuple_to_dict(nt::NamedTuple)
    Dict{String, Float64}([string(p[1])=>Float64(p[2]) for p in pairs(nt)])
end

"""
    fit_spectrum(obs_wls, obs_flux, obs_err, linelist, initial_guesses, fixed_params; kwargs...)

Find the parameters and abundances that best match a rectified observed spectrum.

# Arguments:
- `obs_wls`: the wavelengths of the observed spectrum in Å
- `obs_flux`: the rectified flux of the observed spectrum
- `obs_err`: uncertainty in `flux`
- `linelist`: a linelist to use for the synthesis
- `initial_guesses`: a NamedTuple specifying initial guesses for the parameters to be fit.  See 
  "Specifying parameters" below.
- `fixed_params`: a NamedTuple specifying parameters to be held fixed. See "Specifying parameters" 
  below.

`initial_guesses` and `fixed_params` can also be specified as Dicts instead of NamedTuples, which is 
more convienient when calling Korg from python.

# Specifying parameters
Parameters are specified as NamedTuples, which look like this: `(Teff=5000, logg=4.5, m_H=0.0)`.
Single-element NamedTuples require a semicolon: `(; Teff=5000)`. 
## Required parameters
`Teff` and `logg` *must* be specified in either `initial_guesses` or `fixed_params`.
## Optional Parameters
These can be specified in either `initial_guesses` or `fixed_params`, but if they are not default 
values are used.
- `m_H`: the metallicity of the star, in dex. Default: `0.0`
- `alpha_H`: the alpha enhancement of the star, in dex. Default: `m_H`.  Note that, because of the 
  parameter range supported by [`Korg.interpolate_marcs`](@ref), only values within ±1 of `m_H` 
  are supported.
- `vmic`: the microturbulence velocity, in km/s. Default: `1.0`
- `vsini`: the projected rotational velocity of the star, in km/s. Default: `0.0`. 
   See [`Korg.apply_rotation`](@ref) for details.
- `epsilon`: the linear limb-darkening coefficient. Default: `0.6`. Used for applying rotational 
  broadening only.  See [`Korg.apply_rotation`](@ref) for details.
- Individual elements, e.g. `Na`, specify the solar-relative ([X/H]) abundance of that element. 

!!! tip
    If you are doing more than a few fits, you will save a lot of time by precomputing the LSF 
    matrix and synthesis wavelengths.  See the keyword arguments below for how to do that.

# Keyword arguments
- `windows` (optional) is a vector of wavelength pairs, each of which specifies a wavelength 
  "window" to synthesize and contribute to the total χ². If not specified, the entire spectrum is 
  used. Overlapping windows are automatically merged.
- `LSF_matrix` (optional) is a matrix which maps the synthesized spectrum to the observed spectrum. 
  If not specified, it is calculated using `Korg.compute_LSF_matrix`.  Computing the LSF matrix can 
  be expensive, so you may want to precompute it if you are fitting many spectra with the same LSF.
- `synthesis_wls`: a superset of the wavelengths to synthesize, as a range.  If not specified, 
   wavelengths spanning the first and last windows are used. If you pass in a precomputed 
   LSF matrix, you must make sure that the synthesis wavelengths match it.
- `wl_buffer` is the number of Å to add to each side of the synthesis range for each window.
- `precision` specifies the tolerance for the solver to accept a solution. The solver operates on 
   transformed parameters, so `precision` doesn't translate straitforwardly to Teff, logg, etc, but 
   the default is, `1e-3`, provides a worst-case tolerance of about 1.5K in `Teff`, 0.002 in `logg`, 
   0.001 in `m_H`, and 0.004 in detailed abundances.


# Returns
A NamedTuple with the following fields:
- `best_fit_params`: the best-fit parameters
- `best_fit_flux`: the best-fit flux, with LSF applied, resampled, and rectified.  
- `obs_wl_mask`: a bitmask for `obs_wls` which selects the wavelengths used in the fit (i.e. those 
  in the `windows`)
- `solver_result`: the result object from `Optim.jl`
- `trace`: a vector of NamedTuples, each of which contains the parameters at each step of the 
  optimization.

!!! tip
    The function takes a long time to compile the first time it is called. Compilation performance 
    is significantly better on Julia 1.10 than previous versions, so if you are using an older
    version of Julia, you may want to upgrade.
"""
function fit_spectrum(obs_wls, obs_flux, obs_err, linelist, initial_guesses, fixed_params=(;);
                      windows=[(obs_wls[1], obs_wls[end])],
                      synthesis_wls = obs_wls[1] - 10 : 0.01 : obs_wls[end] + 10,
                      R=nothing, 
                      LSF_matrix = if isnothing(R)
                        throw(ArgumentError("Either R or LSF_matrix must be specified."))
                      else
                          Korg.compute_LSF_matrix(synthesis_wls, obs_wls, R)
                      end,
                      wl_buffer=1.0, precision=1e-3)

    initial_guesses, fixed_params = validate_params(initial_guesses, fixed_params)
    ps = collect(pairs(scale(initial_guesses)))
    params_to_fit = first.(ps)
    p0 = last.(ps) # the initial guess as a vector of scaled values

    @assert length(initial_guesses) > 0 "Must specify at least one parameter to fit."

    # calculate some synth ranges which span all windows, and the LSF submatrix that maps to them only
    windows = merge_bounds(windows, 2wl_buffer)
    obs_wl_mask, synth_wl_mask, multi_synth_wls = 
        calculate_multilocal_masks_and_ranges(windows, obs_wls, synthesis_wls, wl_buffer)

    chi2 = let data=obs_flux[obs_wl_mask], obs_err=obs_err[obs_wl_mask], synthesis_wls=multi_synth_wls, 
               LSF_matrix=LSF_matrix[obs_wl_mask, synth_wl_mask], linelist=linelist, 
               params_to_fit=params_to_fit, fixed_params=fixed_params
        function chi2(scaled_p)
            # this extremely weak prior helps to regularize the optimization
            negative_log_scaled_prior = sum(@. scaled_p^2/100^2)
            guess = unscale(Dict(params_to_fit .=> scaled_p))
            params = merge(guess, fixed_params)
            flux = try
                synthetic_spectrum(synthesis_wls, linelist, LSF_matrix, params)
            catch e
                if e isa Korg.ChemicalEquilibriumError
                    # This is a nice huge chi2 value, but not too big.  It's what you get if 
                    # difference at each pixel in the (rectified) spectra is 1, which is 
                    # more-or-less an upper bound.
                    return sum(1 ./ obs_err.^2) 
                else
                    rethrow(e)
                end
            end
            sum(((flux .- data)./obs_err).^2) + negative_log_scaled_prior
        end
    end 
    
    res = if length(p0) == 1
        # if we are fitting a single parameter, experimentation shows that Nelder-Mead (the default)
        # is faster than BFGS

        # there seems to be a problem with trace storage for this optimizer, so we don't request it
        # the precision keyword is also ignored
        optimize(chi2, p0, Optim.Options(x_tol=precision); autodiff=:forward) 
    else
        # if we are fitting a multiple parameters, use BFGS with autodiff
        optimize(chi2, p0, BFGS(linesearch=LineSearches.BackTracking()),
                 Optim.Options(x_tol=precision, time_limit=10_000, store_trace=true, 
                               extended_trace=true); autodiff=:forward)
    end
    solution = unscale(Dict(params_to_fit .=> res.minimizer))

    trace = map(res.trace) do t
        unscaled_params = unscale(Dict(params_to_fit .=> t.metadata["x"]))
        unscaled_params["chi2"] = t.value
        unscaled_params
    end

    full_solution = merge(solution, fixed_params)
    best_fit_flux = try
        synthetic_spectrum(multi_synth_wls, linelist, LSF_matrix[obs_wl_mask, synth_wl_mask], full_solution)
    catch e
        println(e)
    end

    (best_fit_params=solution, best_fit_flux=best_fit_flux, obs_wl_mask=obs_wl_mask, 
     solver_result=res, trace=trace)
end

"""
Sort a vector of lower-bound, upper-bound pairs and merge overlapping ranges.  Used by 
fit_spectrum.
"""
function merge_bounds(bounds, merge_distance)
    bounds = sort(bounds, by=first)
    new_bounds = [bounds[1]]
    for i in 2:length(bounds)
        # if these bounds are within merge_distance of the previous, extend the previous, 
        # otherwise add them to the list
        if bounds[i][1] <= new_bounds[end][2] + merge_distance 
            new_bounds[end] = (new_bounds[end][1], max(bounds[i][2], new_bounds[end][2]))
        else
            push!(new_bounds, bounds[i])
        end
    end
    new_bounds
end


"""
    calculate_multilocal_masks_and_ranges(obs_bounds_inds, obs_wls, synthesis_wls)

Given a vector of target synthesis ranges in the observbed spectrum, return the masks, etc required.

Arguments:
    - `windows`: a vector of pairs of wavelength lower and upper bounds.
    - `obs_wls`: the wavelengths of the observed spectrum
    - `synthesis_wls`: the wavelengths of the synthesis spectrum
    - `wl_buffer`: the number of Å to add to each side of the synthesis range

Returns:
    - `obs_wl_mask`: a bitmask for `obs_wls` which selects the observed wavelengths
    - `synthesis_wl_mask`: a bitmask for `synthesis_wls` which selects the synthesis wavelengths 
       needed to generated the masked observed spectrum.
    - `multi_synth_wls`: The vector of ranges to pass to `Korg.synthesize`.
"""
function calculate_multilocal_masks_and_ranges(windows, obs_wls, synthesis_wls, wl_buffer)
    # bitmasks for obs_wls synthesis_wls to isolate the subspectra
    obs_wl_mask = zeros(Bool, length(obs_wls)) 
    synth_wl_mask = zeros(Bool, length(synthesis_wls)) 

    # multi_synth_wls is the vector of wavelength ranges that gets passed to synthesize
    multi_synth_wls = map(windows) do (ll, ul)
        lb, ub = (findfirst(obs_wls .>= ll), findlast(obs_wls .<= ul))
        if isnothing(lb) || isnothing(ub) || lb > ub
            error("The range $ll to $ul is not in the observed spectrum")
        end

        obs_wl_mask[lb:ub] .= true

        synth_wl_lb = findfirst(synthesis_wls .>= obs_wls[lb] - wl_buffer)
        synth_wl_ub = findfirst(synthesis_wls .> obs_wls[ub] + wl_buffer) - 1
        synth_wl_mask[synth_wl_lb:synth_wl_ub] .= true

        synthesis_wls[synth_wl_lb:synth_wl_ub]
    end
    obs_wl_mask, synth_wl_mask, multi_synth_wls
end

"""
    linelist_neighbourhood_indices(linelist, ew_window_size)

Group lines together such that no two lines are closer than twice the value of `line_buffer`.

# Arguments:
- `linelist`: A vector of [`Korg.Line`](@ref)s (see [`Korg.read_linelist`](@ref), 
   [`Korg.get_APOGEE_DR17_linelist`](@ref), and [`Korg.get_VALD_solar_linelist`](@ref)).
- `ew_window_size`: the minimum separation (in Å) either side of lines in a group

# Returns
A vector of vectors, where each inner vector contains the indices of lines in a group.
"""
function linelist_neighbourhood_indices(linelist, ew_window_size)
    linelist_neighbourhood_indices = []        
    current_group = [1]    
    ew_window_size_overlap_cm = 2 * 1e-8 * ew_window_size
    for i in 2:length(linelist)
        if (linelist[i].wl - linelist[current_group[end]].wl) > ew_window_size_overlap_cm
            push!(current_group, i)
        else
            push!(linelist_neighbourhood_indices, current_group)
            current_group = [i]  
        end
    end
    push!(linelist_neighbourhood_indices, current_group)
    linelist_neighbourhood_indices
end

"""
    ews_to_abundances(atm, linelist, A_X, measured_EWs; kwargs... )

Compute per-line abundances on the linear part of the curve of growth given a model atmosphere and a
list of lines with equivalent widths.

# Arguments:
- `atm`: the model atmosphere (see [`Korg.read_model_atmosphere`](@ref) and 
  [`Korg.interpolate_marcs`](@ref)).
- `linelist`: A vector of [`Korg.Line`](@ref)s (see [`Korg.read_linelist`](@ref)).  The lines must 
  be sorted by wavelength.
- `A_X`: a vector containing the A(X) abundances (log(n_X/n_H) + 12) for elements from hydrogen to 
  uranium (see [`Korg.format_A_X`](@ref)). All syntheses are done with these abundances, so if the 
  resulting abundances deviate significantly from these, you may wish to iterate.
- `measured_EWs`: a vector of equivalent widths (in mÅ)

# Returns
A vector of abundances (`A(X) = log10(n_X/n_H) + 12` format) for each line in `linelist`.

# Optional arguments:
- `wl_step` (default: 0.01) is the resolution in Å at which to synthesize the spectrum around each 
   line. 
- `ew_window_size` (default: 2): the farthest (in Å) to consider equivalent width contributions for 
   each line.
All other keyword arguments are passed to [`Korg.synthesize`](@ref) when synthesizing each line.
"""
function ews_to_abundances(atm, linelist, A_X, measured_EWs, ew_window_size::Real=2.0, wl_step=0.01; 
                           synthesize_kwargs...)
    synthesize_kwargs = Dict(synthesize_kwargs)
    if get(synthesize_kwargs, :hydrogen_lines, false)
        throw(ArgumentError("hydrogen_lines must be disabled"))
    end

    if length(linelist) != length(measured_EWs)
        throw(ArgumentError("length of linelist does not match length of ews ($(length(linelist)) != $(length(measured_EWs)))"))
    end
    
    if !issorted(linelist; by=l->l.wl) 
        throw(ArgumentError("linelist must be sorted"))
    end

    if any(l -> Korg.ismolecule(l.species), linelist)
        throw(ArgumentError("linelist contains molecular species"))
    end

    # Check that the user is supplying EWs in mA
    if 1 > maximum(measured_EWs)
        @warn "Maximum EW given is less than 1 mA. Check that you're giving EWs in mÅ (*not* Å)."
    end

    # Group lines together ensuring that no λ is closer to it's neighbour than twice the ew_window_size.
    group_indices = linelist_neighbourhood_indices(linelist, ew_window_size)

    A0_minus_log10W0 = Array{Float64}(undef, length(linelist))
    for indices in group_indices
        wl_ranges = map(linelist[indices]) do line
            # constructing the range this way ensures that the last point is exactly λ_stop
            # and there is always a point on the line center
            λ_start, λ_stop = (1e8 * line.wl - ew_window_size, 1e8 * line.wl + ew_window_size)
            range(λ_start, λ_stop; length=Int(round((λ_stop - λ_start)/wl_step))+1)
        end

        spectrum = Korg.synthesize(atm, linelist[indices], A_X, wl_ranges; synthesize_kwargs...)

        for (i, (idx, line)) in enumerate(zip(spectrum.subspectra, linelist[indices]))
            depth = 1 .- spectrum.flux[idx] ./ spectrum.cntm[idx]
            logEW = log10(trapz(spectrum.wavelengths[idx], depth) * 1e3) # convert to mÅ
            A0_minus_log10W0[indices[i]] = A_X[Korg.get_atoms(line.species)[1]] - logEW
        end
    end

    log10.(measured_EWs) .+ A0_minus_log10W0
end

end # module