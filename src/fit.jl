"""
Functions for fitting to data.

!!! warning
    This submodule is in beta. It's API may change.
"""
module Fit
using ..Korg, LineSearches, Optim
using Interpolations: linear_interpolation
using ForwardDiff, DiffResults
using Trapz
using Statistics: mean, std

# used by scale and unscale for some parameters
function tan_scale(p, lower, upper) 
    if !(lower <= p <= upper)
        throw(ArgumentError("p=$p is not in the range $lower to $upper"))
    end
    tan(π * (((p-lower)/(upper-lower))-0.5))
end
tan_unscale(p, lower, upper) = (atan(p)/π + 0.5)*(upper - lower) + lower

# these are the parameters which are scaled by tan_scale
const tan_scale_params = Dict(
    "epsilon" => (0, 1),
    "cntm_offset" => (-0.5, 0.5),
    "cntm_slope" => (-0.1, 0.1),
    # we can't get these directly from Korg.get_atmosphere_archive() because it will fail in the 
    # test environment, but they are simply the boundaries of the SDSS marcs grid used by
    # Korg.interpolate_marcs.
    "Teff" => (2800, 8000),
    "logg" => (-0.5, 5.5),
    "m_H" => (-5, 1),
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
an internal function oned by fitting routines.
See [`Korg.synthesize`](@ref) to synthesize spectra as a Korg user.
"""
function synthetic_spectrum(synthesis_wls, linelist, LSF_matrix, params, synthesis_kwargs)
    specified_abundances = Dict([p for p in pairs(params) if p.first in Korg.atomic_symbols])
    alpha_H = "alpha_H" in keys(params) ? params["alpha_H"] : params["m_H"]
    A_X::Vector{valtype(params)} = Korg.format_A_X(params["m_H"], alpha_H, specified_abundances; solar_relative=true)

    # clamp_abundances clamps M_H, alpha_M, and C_M to be within the atm grid
    atm = Korg.interpolate_marcs(params["Teff"], params["logg"], A_X; clamp_abundances=true, perturb_at_grid_values=true)

    sol = Korg.synthesize(atm, linelist, A_X, synthesis_wls; vmic=params["vmic"], line_buffer=0,
                          electron_number_density_warn_threshold=Inf, synthesis_kwargs...)

    # apply cntm adjustments
    central_wavelength = (sol.wavelengths[begin] + sol.wavelengths[end]) / 2
    cntm_adjustment = 1 .- params["cntm_offset"] .- params["cntm_slope"]*(sol.wavelengths .- central_wavelength)
    F = sol.flux ./ (sol.cntm .* cntm_adjustment)

    F = Korg.apply_rotation(F, synthesis_wls, params["vsini"], params["epsilon"])
    LSF_matrix * F
end

"""
Validate fitting parameters, and insert default values when needed. Used by [`fit_spectrum`](@ref).

these can be specified in either initial_guesses or fixed_params, but if they are not, these values
 are inserted into fixed_params
"""
function validate_params(initial_guesses::AbstractDict, fixed_params::AbstractDict;
                         required_params = ["Teff", "logg"],
                         default_params = Dict("m_H"=>0.0, "vsini"=>0.0, "vmic"=>1.0, "epsilon"=>0.6, "cntm_offset"=>0.0, "cntm_slope"=>0.0),
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
validate_params(initial_guesses::AbstractDict, fixed_params::NamedTuple; kwargs...) = 
    validate_params(initial_guesses, _namedtuple_to_dict(fixed_params); kwargs...)
validate_params(initial_guesses::NamedTuple, fixed_params=AbstractDict{String, Float64}(); kwargs...) = 
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
more convenient when calling Korg from python.

# Specifying parameters
Parameters are specified as named tuples or dictionaries. Named tuples look like this:
 `(Teff=5000, logg=4.5, m_H=0.0)`.  Single-element named tuples require a semicolon: `(; Teff=5000)`.

### Required parameters
`Teff` and `logg` *must* be specified in either `initial_guesses` or `fixed_params`.
### Optional Parameters
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
- `cntm_offset`: a constant offset to the continuum. Default: `0.0` (no correction).
- `cntm_slope`: the slope of a linear correction applied to the continuum. Default: `0.0` (no 
  correction). The linear correction will be zero at the central wavelength. While it's possible to
  fit for a continuum slope with a fixed offset, it is highly disrecommended.

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
   transformed parameters, so `precision` doesn't translate straightforwardly to Teff, logg, etc, but 
   the default value, `1e-4`, provides a theoretical worst-case tolerance of about 0.15 K in `Teff`, 
   0.0002 in `logg`, 0.0001 in `m_H`, and 0.0004 in detailed abundances. In practice the precision 
   achieved by the optimizer is about 10x bigger than this.
Any additional keyword arguments will be passed to [`Korg.synthesize`](@ref) when synthesizing the
spectra for the fit.

# Returns
A NamedTuple with the following fields:
- `best_fit_params`: the best-fit parameters
- `best_fit_flux`: the best-fit flux, with LSF applied, resampled, and rectified.  
- `obs_wl_mask`: a bitmask for `obs_wls` which selects the wavelengths used in the fit (i.e. those 
  in the `windows`)
- `solver_result`: the result object from `Optim.jl`
- `trace`: a vector of NamedTuples, each of which contains the parameters at each step of the 
  optimization. This is empty for single parameter fits, because the underlying solver doesn't 
  supply it.
- `covariance`: a pair `(params, Σ)` where `params` is vector of parameter name (providing an 
  order), and `Σ` is an estimate of the covariance matrix of the parameters.  It is the approximate 
  inverse hessian of the log likelihood at the best-fit parameter calculated by the BGFS algorithm, 
  and should be interpreted with caution. For single-parameter fits (which are done with 
  Nelder-Mead), this is not provided.

!!! tip
    This function takes a long time to compile the first time it is called. Compilation performance 
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
                      wl_buffer=1.0, precision=1e-4, synthesis_kwargs...)
    if length(obs_wls) != length(obs_flux) || length(obs_wls) != length(obs_err)
        throw(ArgumentError("obs_wls, obs_flux, and obs_err must all have the same length."))
    end
    if length(obs_wls) != size(LSF_matrix, 1)
        throw(ArgumentError("the first dimension of LSF_matrix must be the length of obs_wls."))
    end
    if (length(synthesis_wls) != size(LSF_matrix, 2))
        throw(ArgumentError("the second dimension of LSF_matrix must be the length of synthesis_wls " * 
         "If you provided one as a keyword argument, you must also provide the other."))
    end

    initial_guesses, fixed_params = validate_params(initial_guesses, fixed_params)
    ps = collect(pairs(scale(initial_guesses)))
    params_to_fit = first.(ps)
    p0 = last.(ps) # the initial guess as a vector of scaled values

    @assert length(initial_guesses) > 0 "Must specify at least one parameter to fit."

    # calculate some synth ranges which span all windows, and the LSF submatrix that maps to them only
    windows, _ = merge_bounds(windows, 2wl_buffer)
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
                synthetic_spectrum(synthesis_wls, linelist, LSF_matrix, params, synthesis_kwargs)
            catch e
                if (e isa Korg.ChemicalEquilibriumError) || (e isa Korg.LazyMultilinearInterpError)
                    # chemical equilibrium errors happen for a few unphysical model atmospheres

                    # LazyMultilinearInterpError happens when logg is oob for the low-Z atmosphere grid

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
    
    # call optimization library
    res = optimize(chi2, p0, BFGS(linesearch=LineSearches.BackTracking(maxstep=1.0)),
             Optim.Options(x_tol=precision, time_limit=10_000, store_trace=true, 
                           extended_trace=true); autodiff=:forward)

    best_fit_params = unscale(Dict(params_to_fit .=> res.minimizer))

    best_fit_flux = try
        full_solution = merge(best_fit_params, fixed_params)
        synthetic_spectrum(multi_synth_wls, linelist, LSF_matrix[obs_wl_mask, synth_wl_mask], full_solution, synthesis_kwargs)
    catch e
        println(e)
    end

    trace = map(res.trace) do t
        unscaled_params = unscale(Dict(params_to_fit .=> t.metadata["x"]))
        unscaled_params["chi2"] = t.value
        unscaled_params
    end

    invH = let
        # derivate relating the scaled parameters to the unscaled parameters
        # (used to convert the approximate hessian to a covariance matrix in the unscaled params)
        dp_dscaledp = map(res.minimizer, params_to_fit) do scaled_param, param_name
            ForwardDiff.derivative(scaled_param) do scaled_param
                unscale(Dict(param_name=>scaled_param))[param_name]
            end
        end
        # the fact that the scaling is a diagonal operation means that we can do this as an 
        # element-wise product.  If we think of ds/dp as a (diagonal) matrix, this is equivalent to
        # (ds/dp)^T * invH * (ds/dp)
        invH_scaled = res.trace[end].metadata["~inv(H)"]
        invH_scaled .* dp_dscaledp .* dp_dscaledp'
    end

    (best_fit_params=best_fit_params, best_fit_flux=best_fit_flux, obs_wl_mask=obs_wl_mask, 
     solver_result=res, trace=trace, covariance=(params_to_fit, invH))
end

"""
Sort a vector of lower-bound, upper-bound pairs and merge overlapping ranges.  Used by 
fit_spectrum and ews_to_stellar_parameters.

Returns a pair containing:
- a vector of merged bounds
- a vector of vectors of indices of the original bounds which were merged into each merged bound
"""
function merge_bounds(bounds, merge_distance)
    bound_indices = 1:length(bounds)

    # short by lower bound
    s = sortperm(bounds, by=first)
    bounds = bounds[s]
    bound_indices = bound_indices[s]

    new_bounds = [bounds[1]]
    indices = [[bound_indices[1]]]
    for i in 2:length(bounds)
        # if these bounds are within merge_distance of the previous, extend the previous, 
        # otherwise add them to the list
        if bounds[i][1] <= new_bounds[end][2] + merge_distance 
            new_bounds[end] = (new_bounds[end][1], max(bounds[i][2], new_bounds[end][2]))
            push!(indices[end], bound_indices[i])
        else
            push!(new_bounds, bounds[i])
            push!(indices, [bound_indices[i]])
        end
    end
    new_bounds, indices
end


"""
    calculate_multilocal_masks_and_ranges(obs_bounds_inds, obs_wls, synthesis_wls)

Given a vector of target synthesis ranges in the observed spectrum, return the masks, etc required.

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
        synth_wl_ub = findlast(synthesis_wls .<= obs_wls[ub] + wl_buffer)
        synth_wl_mask[synth_wl_lb:synth_wl_ub] .= true

        synthesis_wls[synth_wl_lb:synth_wl_ub]
    end
    obs_wl_mask, synth_wl_mask, multi_synth_wls
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
   each line.  It's very important that this is large enough to include each line entirely.
- `blend_warn_threshold` (default: 0.01) is the minimum absorption between two lines allowed before 
   triggering a warning that they may be blended.
All other keyword arguments are passed to [`Korg.synthesize`](@ref) when synthesizing each line.
"""
function ews_to_abundances(atm, linelist, A_X, measured_EWs; ew_window_size::Real=2.0, wl_step=0.01,
                           blend_warn_threshold=0.01, synthesize_kwargs...)
    synthesize_kwargs = Dict(synthesize_kwargs)
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

    merged_windows, lines_per_window = 
        merge_bounds([(line.wl*1e8 - ew_window_size, line.wl*1e8 + ew_window_size) for line in linelist], 0.0)
    wl_ranges = map(merged_windows) do (wl1, wl2)
        wl1:wl_step:wl2
    end

    # hydrogen_lines should be disabled for most accurate equivalent widths.  This can be overridden
    # by passing hydrogen_lines=true as a keyword argument (included in synthesize_kwargs)
    # line_buffer=0.0 makes things a bit faster, and it causes no problems as long as ew_window_size
    # is sufficient, which is necessary anyway.
    sol = Korg.synthesize(atm, linelist, A_X, wl_ranges; line_buffer=0.0, hydrogen_lines=false, synthesize_kwargs...)
    depth = 1 .- sol.flux ./ sol.cntm

    element_type = promote_type(eltype(A_X), eltype(Korg.get_temps(atm)))
    A0_minus_log10W0 = Array{element_type}(undef, length(linelist))
    all_boundaries = Float64[]
    for (wl_range, subspec, line_indices) in zip(wl_ranges, sol.subspectra, lines_per_window)
        absorption = depth[subspec]

        # get the wl-index of least absorption between each pair of lines
        boundary_indices = map(1:length(line_indices) - 1) do i
            wl1 = linelist[line_indices[i]].wl * 1e8
            wl2 = linelist[line_indices[i+1]].wl * 1e8
            l1_ind = Int(round((wl1 - wl_range[1]) / step(wl_range))) + 1
            l2_ind = Int(round((wl2 - wl_range[1]) / step(wl_range))) + 1
            boundary_index = argmin(absorption[l1_ind:l2_ind]) + l1_ind - 1
            if absorption[boundary_index] > blend_warn_threshold
                @warn "Lines $(line_indices[i]) and $(line_indices[i+1]) ($(linelist[line_indices[i]].wl*1e8) Å and $(linelist[line_indices[i+1]].wl*1e8)) Å appear to be blended.  Between them, the absorption never drops below $(blend_warn_threshold) (minimum: $(ForwardDiff.value(absorption[boundary_index]))). You can adjust this threshold with the blend_warn_threshold keyword argument."
            end
            boundary_index
        end
        boundary_indices = [1 ; boundary_indices ; length(subspec)]
        for b in boundary_indices
            push!(all_boundaries, wl_range[b])
        end

        for i in 1:length(line_indices)
            r = boundary_indices[i]:boundary_indices[i+1]
            logEW = log10(trapz(wl_range[r], absorption[r]) * 1e3) # convert to mÅ
            Z = Korg.get_atoms(linelist[line_indices[i]].species)[1]
            A0_minus_log10W0[line_indices[i]] = A_X[Z] - logEW
        end
    end

    # TODO maybe return this stuff?
    log10.(measured_EWs) .+ A0_minus_log10W0#, (sol.wavelengths, 1 .- depth), all_boundaries
end

"""
    ews_to_stellar_parameters(linelist, measured_EWs, [measured_EW_err]; kwargs...)


Find stellar parameters from equivalent widths the "old fashioned" way.  This function finds the 
values of ``T_\\mathrm{eff}``, ``\\log g``, ``v_{mic}``, and [m/H] which satisfy the following conditions 
(using a Newton-Raphson solver):
- The slope of the abundances of neutral lines with respect to lower excitation potential is zero.
- The difference between the mean abundances of neutral and ionized lines is zero.
- The slope of the abundances of neutral lines with respect to reduced equivalent width is zero.
- The difference between the mean abundances of all lines and the model-atmosphere input [m/H] is zero.
Here the "slope" refers to the slope of a linear fit to the abundances of the lines in question.

# Arguments:
- `linelist`: A vector of [`Korg.Line`](@ref) objects (see [`Korg.read_linelist`](@ref)).  The lines 
   must be sorted by wavelength.
- `measured_EWs`: a vector of equivalent widths (in mÅ).
- `measured_EW_err` (optional): the uncertainty in `measured_EWs`.  If not specified, all lines are 
   assumed to have the same uncertainty. These uncertainties are used when evaluating the equations 
   above, and are propagated to provide uncertainties in the resulting parameters.

# Returns:
A tuple containing:
- the best-fit parameters: `[Teff, logg, vmic, [m/H]]` as a vector
- the statistical uncertainties in the parameters, propagated from the uncertainties in the 
  equivalent widths. This is zero when EW uncertainties are not specified.
- the systematic uncertainties in the parameters, estimated from the scatter the abundances computed 
  from each line not accounted for by the EW uncertainties.

# Keyword arguments:
- `Teff0` (default: 5000.0) is the starting guess for Teff
- `logg0` (default: 3.5) is the starting guess for logg
- `vmic0` (default: 1.0) is the starting guess for vmic. Note that this must be nonzero in order to 
   avoid null derivatives. Very small values are fine.
- `m_H0` (default: 0.0) is the starting guess for [m/H]
- `tolerances` (default: `[1e-3, 1e-3, 1e-3, 1e-3]`) is the tolerance for the residuals each equation
   listed above. The solver stops when all residuals are less than the corresponding tolerance.
- `max_step_sizes` (default: `[1000.0, 1.0, 0.3, 0.5]`) is the maximum step size to take in each 
   parameter direction.  This is used to prevent the solver from taking too large of a step and 
   missing the solution.  Be particularly cautious with the vmic (third) parameter, as the 
   unadjusted step size is often too large.
- `parameter_ranges` (default: `[(2800.0, 8000.0), (-0.5, 5.5), (1e-3, 10.0), (-2.5, 1.0)]`) is the 
   allowed range for each parameter. This is used to prevent the solver from wandering into 
   unphysical parameter space, or outside the range of the MARCS grid supported by 
   [`Korg.interpolate_marcs`](@ref). The default ranges ``T_\\mathrm{eff}``, ``\\log g``,  and [m/H]
   are the widest supported by the MARCS grid. Note that vmic must be nonzero in order to avoid null 
   derivatives.
- `callback`:  is a function which is called at each step of the optimization.  
   It is passed three arguments: 
    - the current values of the parameters
    - the residuals of each equation being solved
    - the abundances of each line computed with the current parameters.
  You can pass a callback function, to e.g. make a plot of the residuals at each step. 
- `max_iterations` (default: 30) is the maximum number of iterations to allow before stopping the 
   optimization.
"""
function ews_to_stellar_parameters(linelist, measured_EWs, measured_EW_err=nothing; 
                                   Teff0=5000.0, logg0=3.5, vmic0=1.0, m_H0=0.0,
                                   tolerances=[1e-3, 1e-3, 1e-4, 1e-3],
                                   max_step_sizes=[1000.0, 1.0, 0.3, 0.5],
                                   parameter_ranges=[extrema.(Korg._sdss_marcs_atmospheres[1][1:2]) 
                                                    ; (1e-3, 10.0) 
                                                    ; (Korg._low_Z_marcs_atmospheres[1][3][1], Korg._sdss_marcs_atmospheres[1][3][end])],
                                   fix_params=[false, false, false, false],
                                   callback=Returns(nothing), max_iterations=30, passed_kwargs...)
    if :vmic in keys(passed_kwargs)
        throw(ArgumentError("vmic must not be specified, because it is a parameter fit by ews_to_stellar_parameters.  Did you mean to specify vmic0, the starting value? See the documentation for ews_to_stellar_parameters if you would like to fix microturbulence to a given value."))
    end
    if length(linelist) != length(measured_EWs) || (!isnothing(measured_EW_err) && (length(linelist) != length(measured_EW_err)))
        throw(ArgumentError("length of linelist does not match length of ews ($(length(linelist)) != $(length(measured_EWs)))"))
    end
    formulas = [line.species.formula for line in linelist]
    if any(Ref(formulas[1]) .!= formulas)
        throw(ArgumentError("All lines must be from the same element."))
    end
    if Korg.ismolecule(linelist[1].species)
        throw(ArgumentError("Cannot do stellar parameter determination with molecular lines."))
    end
    neutrals = [l.species.charge == 0 for l in linelist]
    if (sum(neutrals) < 3) || (sum(.! neutrals) < 1)
        throw(ArgumentError("Must have at least 3 neutral lines and 1 ion line."))
    end
    if vmic0 == 0.0
        throw(ArgumentError("Starting guess for vmic (vmic0) must be nonzero."))
    end
    if any(p[1] >= p[2] for p in parameter_ranges)
        throw(ArgumentError("The lower bound of each parameter must be less than the upper bound."))
    end
    if parameter_ranges[3][1] <= 0.0
        throw(ArgumentError("The lower bound of vmic must be greater than zero. (vmic must be nonzero in order to avoid null derivatives. Very small values are fine.)"))
    end

    # the widest parameter ranges allowed for model atmosphere interp
    atm_lb = first.(Korg._sdss_marcs_atmospheres[1][1:3])
    atm_lb[3] = Korg._low_Z_marcs_atmospheres[1][3][1]
    atm_ub = last.(Korg._sdss_marcs_atmospheres[1][1:3])
    if any(first.(parameter_ranges[[1, 2, 4]]) .< atm_lb) || any(last.(parameter_ranges[[1, 2, 4]]) .> atm_ub)
        throw(ArgumentError("The parameter ranges must be within the range of the MARCS grid"))
    end

    params0 = [Teff0, logg0, vmic0, m_H0]
    params = clamp(params0, first.(parameter_ranges), last.(parameter_ranges))
    for (p, p0, n) in zip(params, params0, ["Teff", "logg", "vmic", "metallicity"])
        if p != p0
            @warn "Initial guess for $n ($p0) has been clamped to $p, to be within the allowed range."
        end
    end

    # set up closure to compute residuals
    get_residuals = (p) -> 
        _stellar_param_equation_residuals(p, linelist, measured_EWs, measured_EW_err, 
                                           fix_params, callback, passed_kwargs)
    iterations = 0
    J_result = DiffResults.JacobianResult(params)
    while true
        J_result = ForwardDiff.jacobian!(J_result, get_residuals, params)
        J = DiffResults.jacobian(J_result)
        residuals = DiffResults.value(J_result)
        if all((abs.(residuals) .< tolerances)[.! fix_params]) # stopping condition
            break
        end
        step = zeros(length(params))
        step[.! fix_params] = - J[.!fix_params, .!fix_params] \ residuals[.!fix_params]
        params += clamp.(step, -max_step_sizes, max_step_sizes)
        params .= clamp.(params, first.(parameter_ranges), last.(parameter_ranges))

        iterations += 1
        if iterations > max_iterations
            @warn "Failed to converge after $max_iterations iterations.  Returning the current guess."
            return params, fill(NaN, 4), fill(NaN, 4)
        end
    end

    # compute uncertainties
    stat_σ_r, sys_σ_r = _stellar_param_residual_uncertainties(params, linelist, measured_EWs, measured_EW_err, passed_kwargs)

    J = DiffResults.jacobian(J_result)[.! fix_params, .! fix_params]
    stat_σ = zeros(4)
    stat_σ[.! fix_params] .= abs.(J \ stat_σ_r[.! fix_params])
    sys_σ = zeros(4)
    sys_σ[.! fix_params] .= abs.(J \ sys_σ_r[.! fix_params])

    params, stat_σ, sys_σ
end

# called by ews_to_stellar_parameters
function _stellar_param_equation_residuals(params, linelist, EW, EW_err, 
                                           fix_params, callback, passed_kwargs)
    A, A_inv_var, neutrals, REWs, Z = 
        _stellar_param_equations_precalculation(params, linelist, EW, EW_err, passed_kwargs)


    teff_residual = _get_slope([line.E_lower for line in linelist[neutrals]],
                               A[neutrals], A_inv_var[neutrals])
    logg_residual = (_weighted_mean(A[neutrals], A_inv_var[neutrals]) -
                     _weighted_mean(A[.! neutrals], A_inv_var[.! neutrals]))
    vmic_residual = _get_slope(REWs, A[neutrals], A_inv_var[neutrals])
    feh_residual = _weighted_mean(A, A_inv_var) - (params[4] + Korg.grevesse_2007_solar_abundances[Z])
    residuals = [teff_residual, logg_residual, vmic_residual, feh_residual]
    residuals .*= .! fix_params # zero out residuals for fixed parameters

    callback(ForwardDiff.value.(params), ForwardDiff.value.(residuals), ForwardDiff.value.(A))
    residuals
end

# called by _stellar_param_equation_residuals
# returns (statistical_uncertainty, systematic_uncertainty)
function _stellar_param_residual_uncertainties(params, linelist, EW, EW_err, passed_kwargs)
    A, A_inv_var, neutrals, REWs, _ = 
        _stellar_param_equations_precalculation(params, linelist, EW, EW_err, passed_kwargs)

    # estimated total (including systematic) err in the abundances of each line
    total_err = std(A)
    total_ivar = ones(length(A)) * total_err^-2

    stat_sigma, total_sigma = map([A_inv_var, total_ivar]) do ivar
        sigma_mean = 1 ./ sqrt(sum(ivar))
        teff_residual_sigma = _get_slope_uncertainty([line.E_lower for line in linelist[neutrals]], ivar[neutrals])
        vmic_residual_sigma = _get_slope_uncertainty(REWs, ivar[neutrals])
        [teff_residual_sigma, sigma_mean, vmic_residual_sigma, sigma_mean]
    end

    if EW_err == ones(length(EW))
        # in the case that EW uncertainties were specified, the residuals were calculated 
        # assuming errs = 1, but we want to ignore them here, and call all error "systematic"
        [0.0, 0.0, 0.0, 0.0], total_sigma
    else
        stat_sigma, sqrt.(max.(total_sigma.^2 .- stat_sigma.^2, 0))
    end
end

function _stellar_param_equations_precalculation(params, linelist, EW, EW_err, passed_kwargs)
    teff, logg, vmic, feh = params
    A_X = Korg.format_A_X(feh)
    atm = Korg.interpolate_marcs(teff, logg, A_X; perturb_at_grid_values=true, clamp_abundances=true)
    A = Korg.Fit.ews_to_abundances(atm, linelist, A_X, EW, vmic=vmic; 
                                                 passed_kwargs...)
    if isnothing(EW_err)
        # no EW uncertainties specified. Let's set the same inverse variance for all lines
        A_inv_var = ones(length(EW))
    else
        # convert error in EW to inverse variance in A (assuming linear part of C.O.G.)
        A_inv_var = (EW .* A ./ EW_err) .^ 2
    end

    neutrals = [l.species.charge == 0 for l in linelist]
    REWs = log10.(EW[neutrals] ./ [line.wl for line in linelist[neutrals]])
    # this is guaranteed not to be a mol (checked by ews_to_stellar_parameters).
    Z = Korg.get_atoms(linelist[1].species)[1]

    A, A_inv_var, neutrals, REWs, Z
end

function _weighted_mean(x, inv_var)
    sum(x .* inv_var) / sum(inv_var)
end

# called by _stellar_param_equation_residuals
function _get_slope(xs, ys, inv_var)
    Δx = xs .- mean(xs)    
    Δy = ys .- mean(ys)
    sum(Δx .* Δy .* inv_var) ./ sum(Δx.^2 .* inv_var)
end

function _get_slope_uncertainty(xs, ivar::AbstractVector) # guard against scalar ivar
    sqrt(sum(ivar) / (sum(ones(length(xs)) .* ivar)*sum(ivar .* xs.^2) - sum(ivar .* xs)^2))
end

end # module
