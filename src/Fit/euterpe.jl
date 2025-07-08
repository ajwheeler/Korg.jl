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
    model_spectra = hcat(model_spectra...) # shape (n_wls, n_models)

    # TODO better
    pseudocontinuum_mask = abs.(model_spectra[:, end] - model_spectra[:, 1]) .< 0.01 # TODO kwarg
    model_pseudocontinua = [calculate_pseudocontinuum(obs_wls, m, pseudocontinuum_mask)
                            for m in eachcol(model_spectra)]
    model_pseudocontinua = hcat(model_pseudocontinua...)
    obs_pseudocontinuum = calculate_pseudocontinuum(obs_wls, obs_flux, pseudocontinuum_mask)

    EWs = map(eachcol(model_spectra), eachcol(model_pseudocontinua)) do m, p
        calculate_rough_EWs(obs_wls, m, p, line_indices)
    end
    EWs = hcat(EWs...)
    obs_EWs = calculate_rough_EWs(obs_wls, obs_flux, obs_pseudocontinuum, line_indices)
    EW_itps = interpolate_and_predict(abundance_perturbations, obs_EWs, EWs)
    EW_As = [itp(EW) for (itp, EW) in zip(EW_itps, obs_EWs)]

    zscores = ((obs_flux .- model_spectra) ./ err) .^ 2
    chi2 = Matrix{Float64}(undef, length(line_indices), length(abundance_perturbations))
    for (i, r) in enumerate(line_indices)
        chi2[i, :] .= sum(zscores[r, :]; dims=1)[:]
    end
    chi2_itps, chi2_coeffs = quadratic_minimizers(abundance_perturbations, chi2)
    chi2_As = -chi2_coeffs[2, :] ./ 2chi2_coeffs[3, :]

    depths = [calculate_line_core_depths(m, pc, line_indices)
              for (m, pc) in zip(eachcol(model_spectra), eachcol(model_pseudocontinua))]
    depths = hcat(depths...)
    obs_depths = calculate_line_core_depths(obs_flux, obs_pseudocontinuum, line_indices)
    depth_itps = interpolate_and_predict(abundance_perturbations, obs_depths, depths)
    depth_As = [itp.(d) for (itp, d) in zip(depth_itps, obs_depths)]

    (; EW_As, depth_As, chi2_As, obs_wl_mask, model_spectra, pseudocontinuum_mask,
     obs_pseudocontinuum, EWs, obs_EWs, EW_itps, line_indices, model_pseudocontinua, chi2,
     chi2_coeffs, chi2_itps, depths, obs_depths, depth_itps)
end

"""
Fit a quadratic model to chi2(abundance_perturbation), and report the minimizer
"""
function quadratic_minimizers(abundance_perturbations, chi2)
    # fit a quadratic model
    A = hcat((abundance_perturbations .^ i for i in 0:2)...)
    coeffs = A \ chi2'
    map(eachcol(coeffs)) do c
        x -> sum(c .* x .^ (0:2))
    end, coeffs
end

function interpolate_and_predict(abundance_perturbations, observed_quantities, model_quantities)
    map(eachrow(model_quantities)) do quantities
        A = hcat((quantities .^ i for i in 0:3)...)
        coeffs = A \ abundance_perturbations
        x -> sum(coeffs .* x .^ (0:3))
    end
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

function calculate_line_core_depths(spectrum, pseudocontinuum, line_indices)
    map(line_indices) do r
        # get the middle 5 pixels in the range
        # TODO everthing about this is terrible
        firstind = r[1] + length(r) ÷ 2 - 2
        lastind = firstind + 4
        firstind = max(1, firstind)
        lastind = min(length(spectrum), lastind)
        mean(spectrum[firstind:lastind] - pseudocontinuum[firstind:lastind])
    end
end

end
