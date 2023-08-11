"""
Functions for fitting to data.

!!! warning
    This submodule is in beta. It's API may change.
"""
module Fit
using ..Korg, LineSearches, Optim
using Interpolations: LinearInterpolation
using ForwardDiff
export find_best_fit_params

# used by scale and unscale for some parameters
function tan_scale(p, lower, upper) 
    if !(lower <= p <= upper)
        raise(ArgumentError("p=$p is not in the range $lower to $upper"))
    end
    tan(π * (((p-lower)/(upper-lower))-0.5))
end
tan_unscale(p, lower, upper) = (atan(p)/π + 0.5)*(upper - lower) + lower

# these are the parmeters which are scaled by tan_scale
const tan_scale_params = Dict(
    :epsilon => (0, 1),
    # we can't get these directly from Korg.get_atmosphere_archive() because it will fail in the 
    # test environent, but they are simply the boundaries of the SDSS marcs grid used by
    # Korg.interpolate_marcs.
    :Teff => (2800, 8000),
    :logg => (-0.5, 5.5),
    :m_H => (-2.5, 1),
    map(Korg.atomic_symbols) do el
        A_X_sun = Korg.default_solar_abundances[Korg.atomic_numbers[el]]
        Symbol(el) => (A_X_sun - 10, A_X_sun + 2)
    end...
)

"""
Rescale each parameter so that it lives on (-∞, ∞).
"""
function scale(params::NamedTuple)
    new_pairs = map(collect(pairs(params))) do (name, p)
        name => if name in keys(tan_scale_params)
            tan_scale(p, tan_scale_params[name]...)
        elseif name in [:vmic, :vsini]
            tan_scale(sqrt(p), 0, sqrt(250))
        else
            @error "$name is not a parameter I know how to scale."
        end
    end
    (; new_pairs...)
end

"""
Unscale each parameter so that it lives on the appropriate range instead of (-∞, ∞).
"""
function unscale(params)
    new_pairs = map(collect(pairs(params))) do (name, p)
        name => if name in keys(tan_scale_params)
            tan_unscale(p, tan_scale_params[name]...)
        elseif name in [:vmic, :vsini]
            tan_unscale(p, 0, sqrt(250))^2
        else
            @error "$name is not a parameter I know how to unscale."
        end
    end
    (; new_pairs...)
end

"""
Synthesize a spectrum, returning the flux, with LSF applied, resampled, and rectified.  This is 
used by fitting routines. See [`Korg.synthesize`](@ref) to synthesize spectra as a Korg user.
"""
function synthetic_spectrum(synthesis_wls, linelist, LSF_matrix, params;
                     line_buffer=10)
    specified_abundances = Dict([String(p.first)=>p.second 
                                                  for p in pairs(params) 
                                                  if String(p.first) in Korg.atomic_symbols])
    alpha_H = :alpha_H in keys(params) ? params.alpha_H : params.m_H
    A_X = Korg.format_A_X(params.m_H, alpha_H, specified_abundances; solar_relative=true)

    # clamp_abundances clamps M_H, alpha_M, and C_M to be within the atm grid
    atm = Korg.interpolate_marcs(params.Teff, params.logg, A_X; clamp_abundances=true, perturb_at_grid_values=true)

    sol = Korg.synthesize(atm, linelist, A_X, synthesis_wls; vmic=params.vmic, line_buffer=line_buffer, 
                        electron_number_density_warn_threshold=1e100)
    F = sol.flux ./ sol.cntm
    F = Korg.apply_rotation(F, synthesis_wls, params.vsini, params.epsilon)
    LSF_matrix * F
end

"""
validate fitting parameters, and insert default values when needed.

these can be specified in either initial_guesses or fixed_params, but if they are not, these values are inserted into fixed_params
"""
function validate_params(initial_guesses, fixed_params;
                         required_params = [:Teff, :logg],
                         default_params = (m_H=0.0, vsini=0.0, vmic=1.0, epsilon=0.6),
                         allowed_params = Set([required_params ; keys(default_params)... ; 
                                               Symbol.(Korg.atomic_symbols)]))
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
    end
    fixed_params = merge((; default_params...), fixed_params) 

    # check that no params are both fixed and initial guesses
    let keys_in_both = collect(keys(initial_guesses) ∩ keys(fixed_params))
        if length(keys_in_both) > 0
            throw(ArgumentError("These parameters: $(keys_in_both) are specified as both initial guesses and fixed params."))
        end 
    end

    initial_guesses, fixed_params
end

"""
    find_best_fit_params(obs_wls, obs_flux, obs_err, linelist, initial_guesses, fixed_params; kwargs...)


!!! warning
    This is in beta. It's API may change.  Use at your own risk.

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
"""
function find_best_fit_params(obs_wls, obs_flux, obs_err, linelist, initial_guesses::NamedTuple, 
                              fixed_params::NamedTuple=(;); windows=[(obs_wls[1], obs_wls[end])],
                              synthesis_wls = obs_wls[1] - 10 : 0.01 : obs_wls[end] + 10,
                              R=nothing, 
                              LSF_matrix = if isnothing(R)
                                throw(ArgumentError("Either R or LSF_matrix must be specified."))
                              else
                                  Korg.compute_LSF_matrix(synthesis_wls, obs_wls, R)
                              end,
                              wl_buffer=1.0, precision=1e-3)

    initial_guesses, fixed_params = validate_params(initial_guesses, fixed_params)
    @assert length(initial_guesses) > 0 "Must specify at least one parameter to fit."

    # calculate some synth ranges which span all windows, and the LSF submatrix that maps to them only
    windows = merge_bounds(windows, 2wl_buffer)
    obs_wl_mask, synth_wl_mask, multi_synth_wls = 
        calculate_multilocal_masks_and_ranges(windows, obs_wls, synthesis_wls, wl_buffer)


    chi2 = let data=obs_flux[obs_wl_mask], obs_err=obs_err[obs_wl_mask], synthesis_wls=multi_synth_wls, 
               LSF_matrix=LSF_matrix[obs_wl_mask, synth_wl_mask], linelist=linelist, fixed_params=fixed_params
        scaled_p -> begin
            negative_log_scaled_prior = sum(@. scaled_p^2/100^2)
            guess = unscale((; zip(keys(initial_guesses), scaled_p)...))
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
    
    # the initial guess must be supplied as a vector to optimize
    p0 = collect(values(scale(initial_guesses)))
    res, solution = if length(p0) == 1
        # if we are fitting a single parameter, experimentation shows that Nelder-Mead (the default)
        # is faster than BFGS
        res = optimize(chi2, p0)
        res, unscale((; first(keys(initial_guesses)) => first(res.minimizer)))
    else
        # if we are fitting a multiple parameters, use BFGS with autodiff
        res = optimize(chi2, p0, BFGS(linesearch=LineSearches.BackTracking()),
                 Optim.Options(x_tol=precision, time_limit=10_000, store_trace=true, 
                   extended_trace=true); autodiff=:forward)
        res, unscale((; zip(keys(initial_guesses), res.minimizer)...))
    end
    full_solution = merge(solution, fixed_params)

    scaled_trace = [(; (keys(initial_guesses) .=> t.metadata["x"])...) for t in res.trace]
    trace = map(scaled_trace, res.trace) do p, t
        (; unscale(p)..., chi2=t.value)
    end

    best_fit_flux = try
        synthetic_spectrum(multi_synth_wls, linelist, LSF_matrix[obs_wl_mask, synth_wl_mask], full_solution)
    catch e
        println(e)
    end

    (best_fit_params=solution, best_fit_flux=best_fit_flux, obs_wl_mask=obs_wl_mask, 
     solver_result=res, trace=trace)
end
# make it possible to use dicts instead of NamedTuples for the python people
find_best_fit_params(obs_wls, obs_flux, obs_err, linelist, initial_guesses::Dict, 
                              fixed_params::Dict=Dict{String, Float64}(); kwargs...) = 
    find_best_fit_params(obs_wls, obs_flux, obs_err, linelist, NamedTuple(pairs(initial_guesses)...), 
                         NamedTuple(pairs(fixed_params)...); kwargs...)

"""
Sort a vector of lower-bound, upper-bound pairs and merge overlapping ranges.  Used by 
find_best_params_multilocally.
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
        if isnothing(lb) || isnothing(ub)
            error("The range $ll to $ul is outside the observed spectrum")
        end

        obs_wl_mask[lb:ub] .= true

        synth_wl_lb = findfirst(synthesis_wls .>= obs_wls[lb] - wl_buffer)
        synth_wl_ub = findfirst(synthesis_wls .> obs_wls[ub] + wl_buffer) - 1
        synth_wl_mask[synth_wl_lb:synth_wl_ub] .= true

        synthesis_wls[synth_wl_lb:synth_wl_ub]
    end
    obs_wl_mask, synth_wl_mask, multi_synth_wls
end

end # module