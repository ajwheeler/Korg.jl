## Computation of Q-factors to provide an estimate of the RV precision for a given spectrum as a function of SNR

## Making the approximation that the spectrum as an ``average'' SNR per pixel when estimating the RV precision
"""
    RV_prec_from_Q(Q,SNR,Npix)

Compute the RV precision from the Q factor, the SNR per pixel, and the number of pixels in the spectrum.

# Arguments
- `Q::Real`: Q factor of the spectrum
- `SNR::Real`: SNR per pixel of the spectrum
- `Npix::Real`: Number of pixels in the spectrum

See also: `Qfactor` and `RV_prec_from_noise`
"""
function RV_prec_from_Q(Q::Real,SNR::Real,Npix::Real)
    return (Korg.c_cgs*1e-2)/(Q*sqrt(Npix)*SNR)
end

## Qfactor from the theoreticsl spectrum and LSF model, assuming uniform efficiency as a function of wavelength
"""
    Qfactor(synth_flux,synth_wl,obs_wl,LSF_mat; msk_lres=nothing)

Compute the Q factor from the high resolution theoretical spectrum, the high resolution wavelength grid, the low resolution (observed) wavelength grid, and the LSF matrix.

Based on work from [Bouchy et al. 2001, A&A, 374, 733](https://ui.adsabs.harvard.edu/abs/2001A%26A...374..733B/abstract).

# Arguments
- `synth_flux::Vector{Real}`: High resolution theoretical spectrum
- `synth_wl`: High resolution wavelength grid
- `obs_wl`: Low resolution wavelength grid
- `LSF_mat`: LSF matrix

# Keyword Arguments
- `msk_lres::Vector{Bool}=nothing`: Mask for the low resolution spectrum to account for pixels masked from the observation

See also: `RV_prec_from_Q` and `RV_prec_from_noise`
"""
function Qfactor(synth_flux,synth_wl,obs_wl,LSF_mat; msk_lres=nothing)
    nvecLSF = dropdims(sum(LSF_mat,dims=2),dims=2) # normalisation for the LSF
    spec_lres = LSF_mat*synth_flux./nvecLSF

    dspec_dlam = zeros(length(synth_flux))
    dspec_dlam[2:end] .= diff(synth_flux)./diff(synth_wl)

    Wvec = ((obs_wl.*(LSF_mat*dspec_dlam)./nvecLSF).^2)./spec_lres
    if isnothing(msk_lres)
        return sqrt(sum(Wvec)/sum(spec_lres))
    else
        return sqrt(sum(Wvec[msk_lres])/sum(spec_lres[msk_lres]))
    end
end

## Theoretical limit on RV precision from the measured noise in the spectrum
"""
    RV_prec_from_noise(spec_hres,wave_hres,wave_lres,LSF_mat,noise_contNorm;msk_lres=nothing)

Compute the RV precision using the true observed uncertainties on the continuum normalized spectrum as a function of wavelength.

# Arguments
- `spec_hres::Vector{Real}`: High resolution theoretical spectrum
- `wave_hres`: High resolution wavelength grid
- `wave_lres`: Low resolution wavelength grid
- `LSF_mat`: LSF matrix
- `noise_contNorm::Vector{Real}`: Noise in the continuum normalized spectrum

# Keyword Arguments
- `msk_lres=nothing`: Mask for the low resolution spectrum to account for pixels masked from the observation

See also: `RV_prec_from_Q` and `Qfactor`
"""
function RV_prec_from_noise(spec_hres,wave_hres,wave_lres,LSF_mat,noise_contNorm; msk_lres=nothing)
    nvecLSF = dropdims(sum(LSF_mat,dims=2),dims=2)

    dspec_dlam = zeros(length(spec_hres))
    dspec_dlam[2:end] .= diff(spec_hres)./diff(wave_hres)

    Wvec = ((wave_lres.*(LSF_mat*dspec_dlam)./nvecLSF).^2)./(noise_contNorm.^2)
    if isnothing(msk_lres)
        return (Korg.c_cgs*1e-2) /sqrt(sum(Wvec))
    else
        return (Korg.c_cgs*1e-2) /sqrt(sum(Wvec[msk_lres]))
    end
end