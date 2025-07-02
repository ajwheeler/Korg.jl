"""
Kinda like Bacchus.

TODO
"""
module Euterpe
using ..Korg
using Statistics: std
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
    obs_wls = obs_wls[obs_wl_mask]

    model_spectra = map(abundance_perturbations) do ΔA
        # TODO symbol/str nonsense
        element_specification = Dict(Symbol(element) => ΔA) # TODO solar is implicit here
        @time LSF * synth(; wavelengths=synthesis_wls, element_specification..., synth_kwargs...)[2]
    end
    model_spectra = hcat(model_spectra...)

    # TODO better
    pseudocontinuum_mask = (model_spectra[:, end] - model_spectra[:, 1]) .< 0.01 # TODO kwarg

    @time rough_EWs = map(eachcol(model_spectra)) do model_spectrum
        calculate_rough_EWs(obs_wls, model_spectrum, pseudocontinuum_mask, lines_to_fit)
    end
    @time obs_rough_EWs = calculate_rough_EWs(obs_wls, obs_flux, pseudocontinuum_mask, lines_to_fit)

    (; obs_wl_mask, model_spectra, pseudocontinuum_mask, rough_EWs, obs_rough_EWs)
end

function calculate_rough_EWs(obs_wls, spectrum, pseudocontinuum_mask, lines_to_fit)
    pseudocontinuum = linear_interpolation(obs_wls[pseudocontinuum_mask],
                                           spectrum[pseudocontinuum_mask]; extrapolation_bc=1.0)(obs_wls)
    absorption = pseudocontinuum - spectrum
    map(lines_to_fit) do λ
        # TODO THIS IS DUMB
        r = searchsortedfirst(obs_wls, λ - 0.2):searchsortedfirst(obs_wls, λ + 0.2)
        trapz(obs_wls[r], absorption[r]) * 1e3 # convert to mÅ
    end
end

end
