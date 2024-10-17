function synth(;
               Teff=5000,
               logg=4.5,
               metals_H=0.0,
               alpha_H=0.0,
               abundances=Dict(),
               linelist=get_VALD_solar_linelist(),
               wavelengths=(5000, 6000),
               rectify=true,
               R=Inf,
               vsini=0,
               synthesize_kwargs...)
    A_X = format_A_X(metals_H, alpha_H, abundances)
    atm = interpolate_marcs(Teff, logg, A_X)
    spectrum = synthesize(atm, linelist, A_X, wavelengths; synthesize_kwargs...)
    flux = if rectify
        spectrum.flux ./ spectrum.cntm
    else
        spectrum.flux
    end
    if isfinite(R)
        flux = apply_LSF(flux, spectrum.wavelengths, R)
    end
    if vsini > 0
        flux = apply_rotation(flux, spectrum.wavelengths, vsini)
    end
    spectrum.wavelengths, flux
end