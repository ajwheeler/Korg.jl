## Computation of Q-factors to provide an estimate of the RV precision for a given spectrum as a function of SNR

## Making the approximation that the spectrum as an ``average'' SNR per pixel when estimating the RV precision
function precRV_fromQ(Q,SNR,Npix;c = 299792458) # c in m/s, SNR per pixel
    return c/(Q*sqrt(Npix)*SNR)
end

## Qfactor from the theoreticsl spectrum and LSF model, assuming uniform efficiency as a function of wavelength
function Qfactor(spec_hres,wave_hres,wave_lres,LSF_mat; msk_lres=nothing)
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
function precRV_MeasuredNoise(spec_hres,wave_hres,wave_lres,LSF_mat,noise_contNorm;msk_lres=nothing,c = 299792458)
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