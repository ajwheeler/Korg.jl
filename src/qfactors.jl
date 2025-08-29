"""
    RV_prec_from_Q(Q, RMS_SNR, Npix)

Compute the RV precision (in m/s) from the Q factor, the RMS SNR per pixel, and the number of pixels
in the spectrum.

# Arguments

  - `Q::Real`: Q factor of the spectrum
  - `RMS_SNR::Real`: The root mean squared per pixel SNR of the spectrum. For spectra which are not
    line-blanketed, this is approximately equal to the average SNR.
  - `Npix::Real`: Number of pixels in the spectrum

See also: [`Qfactor`](@ref) and [`RV_prec_from_noise`](@ref)
"""
function RV_prec_from_Q(Q::Real, RMS_SNR::Real, Npix::Real)
    return (Korg.c_cgs * 1e-2) / (Q * sqrt(Npix) * RMS_SNR)
end

"""
    Qfactor(synth_flux, synth_wl, obs_wl, LSF_mat; obs_mask=nothing)

Compute the Q factor from the high-resolution theoretical spectrum, the high-resolution wavelength
grid, the low-resolution (observed) wavelength grid, and the LSF matrix.
Based on work from [Bouchy et al. 2001, A&A, 374, 733](https://ui.adsabs.harvard.edu/abs/2001A%26A...374..733B/abstract).
Note that the Q factor is an approximation when the flux uncertainty is not photon-dominated.

# Arguments

  - `synth_flux`: High-resolution theoretical spectrum
  - `synth_wl`: High-resolution wavelength grid in [any format supported by Korg](@ref wldocs)
  - `obs_wl`: Low-resolution wavelength grid
  - `LSF_mat`: LSF matrix (see [`compute_LSF_matrix`](@ref))

# Keyword Arguments

  - `obs_mask::Vector{Bool}=nothing`: Mask for the low-resolution spectrum to account for pixels masked from the observation

See also: [`RV_prec_from_Q`](@ref) and [`RV_prec_from_noise`](@ref)
"""
function Qfactor(synth_flux, synth_wl, obs_wl, LSF_mat; obs_mask=nothing)
    synth_wl = Wavelengths(synth_wl)
    nvecLSF = dropdims(sum(LSF_mat; dims=2); dims=2) # normalisation for the LSF
    spec_lres = LSF_mat * synth_flux ./ nvecLSF

    dspec_dlam = zeros(length(synth_flux))
    dspec_dlam[2:end] .= diff(synth_flux) ./ diff(synth_wl) * 1e-8 # cm to Å

    Wvec = ((obs_wl .* (LSF_mat * dspec_dlam) ./ nvecLSF) .^ 2) ./ spec_lres
    if isnothing(obs_mask)
        return sqrt(sum(Wvec) / sum(spec_lres))
    else
        return sqrt(sum(Wvec[obs_mask]) / sum(spec_lres[obs_mask]))
    end
end

"""
    RV_prec_from_noise(synth_flux, synth_wl, obs_wl, LSF_mat, obs_err; obs_mask=nothing)

Compute the best achievable RV precision given a spectrum with uncertainties.

# Arguments

  - `synth_flux`: High-resolution theoretical spectrum
  - `synth_wl`: High-resolution wavelength grid in [any format supported by Korg](@ref wldocs)
  - `obs_wl`: Low-resolution wavelength grid
  - `LSF_mat`: LSF matrix (see [`compute_LSF_matrix`](@ref))
  - `obs_err`: Noise in the continuum-normalized spectrum

# Keyword Arguments

  - `obs_mask=nothing`: Mask for the low-resolution spectrum to account for pixels masked from the observation

See also: [`RV_prec_from_Q`](@ref) and [`Qfactor`](@ref)
"""
function RV_prec_from_noise(synth_flux, synth_wl, obs_wl, LSF_mat, obs_err; obs_mask=nothing)
    synth_wl = Wavelengths(synth_wl)
    nvecLSF = dropdims(sum(LSF_mat; dims=2); dims=2)

    dspec_dlam = zeros(length(synth_flux))
    dspec_dlam[2:end] .= diff(synth_flux) ./ diff(synth_wl) * 1e-8 # cm to Å

    Wvec = ((obs_wl .* (LSF_mat * dspec_dlam) ./ nvecLSF) .^ 2) ./ (obs_err .^ 2)
    if isnothing(obs_mask)
        return (Korg.c_cgs * 1e-2) / sqrt(sum(Wvec))
    else
        return (Korg.c_cgs * 1e-2) / sqrt(sum(Wvec[obs_mask]))
    end
end
