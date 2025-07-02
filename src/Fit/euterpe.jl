"""
Kinda like Bacchus.

TODO
"""
module Euterpe
using ..Korg
using Statistics: mean, std
using Interpolations: linear_interpolation
using Trapz: trapz

"""
TODO
TODO don't hardcode Fe / [m/H]
TODO don't hardcode solar abundances
"""
function abund(linelist, element, lines_to_fit, obs_flux, obs_wls, R,
               err=ones(length(obs_flux));
               abundance_perturbations=-0.6:0.3:0.6,
               synth_kwargs...)
    # TODO check no conflicting R specification

    # wavelengths, windows and LSF, TODO refactor
    windows = [(l - 1, l + 1) for l in lines_to_fit] # TODO kwarg
    merged_windows, lines_per_window = Korg.merge_bounds(windows, 0.0)
    synthesis_wls, obs_wl_mask, LSF = Korg.Fit._setup_wavelengths_and_LSF(obs_wls,
                                                                          nothing,
                                                                          nothing, R,
                                                                          merged_windows,
                                                                          0)
    obs_flux = obs_flux[obs_wl_mask]
    err = err[obs_wl_mask]
    obs_wls = obs_wls[obs_wl_mask]
    line_indices = map(lines_to_fit) do λ
        # TODO this is dumb
        searchsortedfirst(obs_wls, λ - 0.2):searchsortedfirst(obs_wls, λ + 0.2)
    end

    model_spectra = map(abundance_perturbations) do ΔA
        # TODO symbol/str nonsense
        element_specification = Dict(Symbol(element) => ΔA) # TODO solar is implicit here
        @time LSF * synth(; wavelengths=synthesis_wls, element_specification..., synth_kwargs...)[2]
    end

    # TODO better
    pseudocontinuum_mask = abs.(model_spectra[end] - model_spectra[1]) .< 0.03 # TODO kwarg
    model_pseudocontinua = calculate_pseudocontinuum.(Ref(obs_wls), model_spectra,
                                                      Ref(pseudocontinuum_mask))
    obs_pseudocontinuum = calculate_pseudocontinuum(obs_wls, obs_flux, pseudocontinuum_mask)

    rough_EWs = calculate_rough_EWs.(Ref(obs_wls), model_spectra, model_pseudocontinua,
                                     Ref(line_indices))
    obs_rough_EW = calculate_rough_EWs(obs_wls, obs_flux, obs_pseudocontinuum, line_indices)

    chi2 = [sum(((obs_flux .- m) / err) .^ 2) for m in model_spectra]

    model_depths = calculate_line_core_depths.(model_spectra, Ref(line_indices))
    obs_depths = calculate_line_core_depths(obs_flux, line_indices)

    (; obs_wl_mask, model_spectra, pseudocontinuum_mask, rough_EWs, obs_rough_EW, line_indices,
     model_pseudocontinua, chi2, model_depths, obs_depths)
end

function calculate_pseudocontinuum(obs_wls, flux, pseudocontinuum_mask)
    linear_interpolation(obs_wls[pseudocontinuum_mask], flux[pseudocontinuum_mask];
                         extrapolation_bc=1.0)(obs_wls)
end

function calculate_rough_EWs(obs_wls, spectrum, pseudocontinuum, line_indices)
    absorption = pseudocontinuum - spectrum
    map(line_indices) do r
        trapz(obs_wls[r], absorption[r]) * 1e3 # convert to mÅ
    end
end

function calculate_line_core_depths(spectrum, line_indices)
    map(line_indices) do r
        # get the middle 5 pixels in the range
        # TODO everthing about this is terrible
        firstind = r[1] + length(r) ÷ 2 - 2
        lastind = firstind + 4
        firstind = max(1, firstind)
        lastind = min(length(spectrum), lastind)
        mean(spectrum[firstind:lastind])
    end
end

end
