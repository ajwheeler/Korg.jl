"""
Functions for fitting to data.

!!! warning
    This submodule is in beta. It's API may change.
"""
module Fit
using ..Korg, LineSearches, Optim
using Interpolations: LinearInterpolation

# TODO specify versions in Project.toml
# TODO derivatives are zero at grid points (perturbing in one direction would be bad???)

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
    map(enumerate([:Teff, :logg, :m_H])) do (ind, p)
        lower = first(Korg.get_atmosphere_archive()[1][ind])
        upper = last(Korg.get_atmosphere_archive()[1][ind])
        p => (lower, upper)
    end...,
    map(Korg.atomic_symbols) do el
        A_X_sun = Korg.default_solar_abundances[Korg.atomic_numbers[el]]
        Symbol(el) => (A_X_sun - 10, A_X_sun + 2)
    end...
)

"""
Rescale each parameter so that it lives on (-∞, ∞) instead of the finite range of the atmosphere 
grid.
"""
function scale(params::NamedTuple)
    new_pairs = map(collect(pairs(params))) do (name, p)
        name => if name in keys(tan_scale_params)
            tan_scale(p, tan_scale_params[name]...)
        elseif name in [:vmic, :vsini]
            sqrt(p)
        else
            @error "$name is not a parameter I know how to scale."
        end
    end
    (; new_pairs...)
end

"""
Unscale each parameter so that it lives on the finite range of the atmosphere grid instead of 
(-∞, ∞).
"""
function unscale(params)
    new_pairs = map(collect(pairs(params))) do (name, p)
        name => if name in keys(tan_scale_params)
            tan_unscale(p, tan_scale_params[name]...)
        elseif name in [:vmic, :vsini]
            p^2
        else
            @error "$name is not a parameter I know how to unscale."
        end
    end
    (; new_pairs...)
end

"""
Synthesize a spectrum, returning the flux, with LSF applied, resampled, and rectified.  This is 
used by fitting routines. See [`Korg.synthesize`](@ref) to synthesize spectra as a Korg user.

TODO remove obs_wls ?
"""
function synthetic_spectrum(synthesis_wls, linelist, LSF_matrix, params;
                     line_buffer=10) :: Vector{<:Real}
    # interpolate_marcs (really deeper in the atm code) wants the types of its args to be the same
    # TODO maybe that requirement should just be lifted instead.

    specified_abundances = Dict([String(p.first)=>p.second for p in pairs(params) if String(p.first) in Korg.atomic_symbols])
    alpha_H = :alpha_H in keys(params) ? params.alpha_H : params.m_H
    A_X = Korg.format_A_X(params.m_H, alpha_H, specified_abundances; solar_relative=false)

    # clamp_abundances clamps M_H, alpha_M, and C_M to be within the atm grid
    atm = Korg.interpolate_marcs(params.Teff, params.logg, A_X; clamp_abundances=true)

    F = try
        Korg.synthesize(atm, linelist, A_X, synthesis_wls; vmic=params.vmic, line_buffer=line_buffer).flux
    catch e
        ones(sum(length.(synthesis_wls)))
    end
    if params.vsini != 0
        F .= apply_rotation(F, synthesis_wls, params.vsini)
    end
    LSF_matrix * F
end

"""
validate fitting parameters, and insert default values when needed.
"""
function validate_fitting_params(initial_guesses, fixed_params)
    required_params = [:Teff, :logg] # these must be specified in either initial_guesses or fixed_params
    let all_params = keys(initial_guesses) ∪ keys(fixed_params)
        for param in required_params
            if !(param in all_params)
                throw(ArgumentError("Must specify $param in either starting_params or fixed_params. (Did you get the capitalization right?)"))
            end
        end
    end

    # these can be specified in either initial_guesses or fixed_params, but if they are not, 
    # these values are inserted into fixed_params
    default_params = (m_H=0.0, vsini=0.0, vmic=1.0)
    default_params = filter(collect(pairs(default_params))) do (k, v)
        !(k in keys(initial_guesses)) && !(k in keys(fixed_params))
    end
    fixed_params = merge((; default_params...), fixed_params) 

    # check that no params are both fixed and initial guesses
    let keys_in_both = collect(keys(initial_guesses) ∩ keys(fixed_params))
        if length(keys_in_both) > 0
            throw(ArgumentError("Theses parameters: $(keys_in_both) are specified as both initial guesses and fixed params."))
        end 
    end

    # TODO check for unknown keys?

    initial_guesses, fixed_params
end

"""
TODO

Given a list of windows in the form of wavelength pairs, find the best stellar parameters by 
synthesizing and fitting within them.

TODO what is buffer and how should it be set?
TODO fit or fix limb darkening epsilon 
TODO test case where windows is unspecified
"""
function find_best_params_windowed(obs_wls, obs_flux, obs_err, linelist, synthesis_wls, 
                                   LSF_matrix, initial_guesses, fixed_params=(;);
                                   windows=[(obs_wls[1], obs_wls[end])], buffer=1.0)
    initial_guesses, fixed_params = validate_fitting_params(initial_guesses, fixed_params)
    @assert length(initial_guesses) > 0 "Must specify at least one parameter to fit."

    # calculate some synth ranges which span all windows, and the LSF submatrix that maps to them only
    windows = merge_bounds(windows, 2buffer)
    obs_wl_inds = map(windows) do (ll, ul)
        (findfirst(obs_wls .> ll), findfirst(obs_wls .> ul)-1)
    end
    obs_wl_mask, masked_obs_wls_range_inds, synth_wl_mask, multi_synth_wls = 
        calculate_multilocal_masks_and_ranges(obs_wl_inds, obs_wls, synthesis_wls, buffer)

    chi2 = let obs_wls=obs_wls[obs_wl_mask], data=obs_flux[obs_wl_mask], 
               obs_err=obs_err[obs_wl_mask], synthesis_wls=multi_synth_wls, 
               LSF_matrix=LSF_matrix[obs_wl_mask, synth_wl_mask], linelist=linelist,
               fixed_params=fixed_params
        scaled_p -> begin
            guess = unscale((; zip(keys(initial_guesses), scaled_p)...))
            params = merge(guess, fixed_params)
            flux = synthetic_spectrum(synthesis_wls, linelist, LSF_matrix, params)
            cntm = synthetic_spectrum(synthesis_wls, [], LSF_matrix, params)
            flux ./= cntm
            sum(((flux .- data)./obs_err).^2)
        end
    end 
    
    # the initial guess must be supplied as a vector to optimize
    p0 = collect(values(scale(initial_guesses)))
    res = optimize(chi2, p0, BFGS(linesearch=LineSearches.BackTracking()),  
                   Optim.Options(x_tol=1e-5, time_limit=10_000, store_trace=true, 
                   extended_trace=true), autodiff=:forward)

    solution = unscale((; zip(keys(initial_guesses), res.minimizer)...))
    full_solution = merge(solution, fixed_params)
    flux = synthetic_spectrum(multi_synth_wls,  linelist,
                 LSF_matrix[obs_wl_mask, synth_wl_mask], full_solution)
    cntm = synthetic_spectrum(multi_synth_wls, [],
                 LSF_matrix[obs_wl_mask, synth_wl_mask], full_solution)
    best_fit_flux = flux ./ cntm

    solution, res, obs_wl_mask, best_fit_flux
end

function find_best_fit_abundance(obs_wls, obs_flux, obs_err, linelist, synthesis_wls, LSF_matrix, 
                                windows, element, initial_guess,fixed_params
                                ; buffer=1.0)
    obs_flux, obs_err, obs_wls, LSF_matrix = mask_out_nans(obs_flux, obs_err, obs_wls, LSF_matrix)
    _, fixed_params = validate_fitting_params((;), fixed_params)

    # TODO if this code remains duplicated between fitting functions, this shoudl be split into
    # a method or modification of calculate_multilocal_masks_and_ranges.  (the merging probably 
    # should as well.)
    # calculate some synth ranges which span all windows, and the LSF submatrix that maps to them only
    windows = merge_bounds(windows, 2buffer)
    obs_wl_inds = map(windows) do (ll, ul)
        (findfirst(obs_wls .> ll), findfirst(obs_wls .> ul)-1)
    end
    obs_wl_mask, masked_obs_wls_range_inds, synth_wl_mask, multi_synth_wls = 
        calculate_multilocal_masks_and_ranges(obs_wl_inds, obs_wls, synthesis_wls, buffer)

    chi2 = let obs_wls=obs_wls[obs_wl_mask], data=obs_flux[obs_wl_mask], 
               obs_err=obs_err[obs_wl_mask], synthesis_wls=multi_synth_wls, 
               LSF_matrix=LSF_matrix[obs_wl_mask, synth_wl_mask], linelist=linelist,
               element=element, fixed_params=fixed_params
               
        scaled_abund -> begin
            abund = unscale((; Symbol(element) => scaled_abund[1]))
            flux = synthetic_spectrum(synthesis_wls, linelist, LSF_matrix, merge(abund, fixed_params))
            cntm = synthetic_spectrum(synthesis_wls, [], LSF_matrix, merge(abund, fixed_params))
            flux ./= cntm
            sum(((flux .- data)./obs_err).^2)
        end
    end
    
    # put scaled initial abundance guess in singleton vector
    p0 = collect(scale((; Symbol(element)=>initial_guess)))
    res = optimize(chi2, p0, autodiff=:forward,
                   Optim.Options(time_limit=10_000, store_trace=true, extended_trace=true))

    best_fit_abund = unscale((; Symbol(element)=>res.minimizer[1]))
    flux = synthetic_spectrum(synthesis_wls, linelist, LSF_matrix, merge(best_fit_abund, fixed_params))
    cntm = synthetic_spectrum(synthesis_wls, [], LSF_matrix, merge(best_fit_abund, fixed_params))
    best_fit_flux = flux ./ cntm

    best_fit_abund[1], res, obs_wl_mask, best_fit_flux
end

"""
Sort a vector of lower-bound, uppoer-bound pairs and merge overlapping ranges.  Used by 
find_best_params_multilocally.
"""
function merge_bounds(bounds, merge_distance)
    bounds = sort(bounds, by=first)
    new_bounds = [bounds[1]]
    for i in 2:length(bounds)
        # if these bounds are within merge_distance of the previous, extend the previous, 
        # otherwise add them to the list
        if bounds[i][1] < new_bounds[end][2] + merge_distance 
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
    - `obs_bounds_inds`: a vector of tuples of the form `(lower_bound_index, upper_bound_index)`
    - `obs_wls`: the wavelengths of the observed spectrum
    - `synthesis_wls`: the wavelengths of the synthesis spectrum
    - `wl_buffer`: the number of Å to add to each side of the synthesis range

Returns:
    - `obs_wl_mask`: a bitmask for `obs_wls` which selects the observed wavelengths
    - `masked_obs_wls_range_inds`: a vector of ranges which select each subspectrum in the masked 
       observed spectrum. 
    - `synthesis_wl_mask`: a bitmask for `synthesis_wls` which selects the synthesis wavelengths 
       needed to generated the masked observed spectrum.
    - `multi_synth_wls`: The vector of ranges to pass to `Korg.synthesize`.
"""
function calculate_multilocal_masks_and_ranges(obs_bounds_inds, obs_wls, synthesis_wls, wl_buffer)
    obs_wl_mask = zeros(Bool, length(obs_wls)) 
    for (lb, ub) in obs_bounds_inds
        obs_wl_mask[lb:ub] .= true
    end
    # these ranges each select a subsprectrum from obs_wls[obs_wl_mask]
    masked_obs_wls_range_inds = [(1 : obs_bounds_inds[1][2] - obs_bounds_inds[1][1] + 1)]
    for (lb, ub) in obs_bounds_inds[2:end]
        prev_ub = masked_obs_wls_range_inds[end][end]
        push!(masked_obs_wls_range_inds, prev_ub+1 : prev_ub+ub-lb+1)
    end
    # bitmask for synthesis_wls to isolate the subspectra
    synth_wl_mask = zeros(Bool, length(synthesis_wls)) 
    # multi_synth_wls is the vector of wavelength ranges that gets passed to synthesize
    multi_synth_wls = map(obs_bounds_inds) do (lb, ub)
        synth_wl_lb = findfirst(synthesis_wls .> obs_wls[lb] - wl_buffer)
        synth_wl_ub = findfirst(synthesis_wls .> obs_wls[ub] + wl_buffer) - 1
        synth_wl_mask[synth_wl_lb:synth_wl_ub] .= true
        synthesis_wls[synth_wl_lb:synth_wl_ub]
    end
    obs_wl_mask, masked_obs_wls_range_inds, synth_wl_mask, multi_synth_wls
end

"""
TODO
"""
function mask_out_nans(obs_flux, obs_err, obs_wls, LSF_matrix)
    mask = @. !isnan(obs_flux) && !isnan(obs_err)
    obs_flux[mask], obs_err[mask], obs_wls[mask], LSF_matrix[mask, :]
end

end # module
