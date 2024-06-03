## Computation of Q-factors to provide an estimate of the RV precision for a given spectrum as a function of SNR

## Making the approximation that the spectrum as an ``average'' SNR per pixel when estimating the RV precision
"""
    precRV_fromQ(Q,SNR,Npix;c = 299792458)

Compute the RV precision from the Q factor, the SNR per pixel, and the number of pixels in the spectrum.

# Arguments
- `Q::Real`: Q factor of the spectrum
- `SNR::Real`: SNR per pixel of the spectrum
- `Npix::Real`: Number of pixels in the spectrum

# Keyword Arguments
- `c::Real=299792458`: Speed of light in m/s
"""
function precRV_fromQ(Q::Real,SNR::Real,Npix::Real;c::Real = 299792458) # c in m/s, SNR per pixel
    return c/(Q*sqrt(Npix)*SNR)
end

## Qfactor from the theoreticsl spectrum and LSF model, assuming uniform efficiency as a function of wavelength
"""
    Qfactor(spec_hres,wave_hres,wave_lres,LSF_mat; msk_lres=nothing)

Compute the Q factor from the high resolution theoretical spectrum, the high resolution wavelength grid, the low resolution (observed) wavelength grid, and the LSF matrix.

# Arguments
- `spec_hres::Vector{Real}`: High resolution theoretical spectrum
- `wave_hres::Vector{Real}`: High resolution wavelength grid
- `wave_lres::Vector{Real}`: Low resolution wavelength grid
- `LSF_mat::Matrix{Real}`: LSF matrix

# Keyword Arguments
- `msk_lres::Vector{Bool}=nothing`: Mask for the low resolution spectrum to account for pixels masked from the observation
"""
function Qfactor(spec_hres::Vector{Real},wave_hres::Vector{Real},wave_lres::Vector{Real},LSF_mat::Matrix{Real}; msk_lres::Vector{Bool}=nothing)
    nvecLSF = dropdims(sum(LSF_mat,dims=2),dims=2) # normalisation for the LSF
    spec_lres = LSF_mat*spec_hres./nvecLSF

    dspec_dlam = zeros(length(spec_hres))
    dspec_dlam[2:end] .= diff(spec_hres)./diff(wave_hres)

    Wvec = ((wave_lres.*(LSF_mat*dspec_dlam)./nvecLSF).^2)./spec_lres
    if isnothing(msk_lres)
        return sqrt(sum(Wvec)/sum(spec_lres))
    else
        return sqrt(sum(Wvec[msk_lres])/sum(spec_lres[msk_lres]))
    end
end

## Theoretical limit on RV precision from the measured noise in the spectrum
"""
    precRV_MeasuredNoise(spec_hres,wave_hres,wave_lres,LSF_mat,noise_contNorm;msk_lres=nothing,c = 299792458)

Compute the RV precision using the true observed uncertainties on the continuum normalized spectrum as a function of wavelength.

# Arguments
- `spec_hres::Vector{Real}`: High resolution theoretical spectrum
- `wave_hres::Vector{Real}`: High resolution wavelength grid
- `wave_lres::Vector{Real}`: Low resolution wavelength grid
- `LSF_mat::Matrix{Real}`: LSF matrix
- `noise_contNorm::Vector{Real}`: Noise in the continuum normalized spectrum

# Keyword Arguments
- `msk_lres::Vector{Bool}=nothing`: Mask for the low resolution spectrum to account for pixels masked from the observation
- `c::Real=299792458`: Speed of light in m/s
"""
function precRV_MeasuredNoise(spec_hres::Vector{Real},wave_hres::Vector{Real},wave_lres::Vector{Real},LSF_mat::Matrix{Real},noise_contNorm::Vector{Real};msk_lres::Vector{Bool}=nothing,c::Real = 299792458)
    nvecLSF = dropdims(sum(LSF_mat,dims=2),dims=2)

    dspec_dlam = zeros(length(spec_hres))
    dspec_dlam[2:end] .= diff(spec_hres)./diff(wave_hres)

    Wvec = ((wave_lres.*(LSF_mat*dspec_dlam)./nvecLSF).^2)./(noise_contNorm.^2)
    if isnothing(msk_lres)
        return c /sqrt(sum(Wvec))
    else
        return c /sqrt(sum(Wvec[msk_lres]))
    end
end