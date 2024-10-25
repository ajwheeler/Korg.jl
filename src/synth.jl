function synth(;
               Teff=5000,
               logg=4.5,
               metals_H=0.0,
               alpha_H=0.0,
               linelist=get_VALD_solar_linelist(),
               wavelengths=(5000, 6000),
               rectify=true,
               R=Inf,
               vsini=0,
               kwargs...,)
    # if a keyword argument is an atomic symbol, it is assumed to be an abundance
    # otherwise, it is passed to synthesize
    abunds_set = [k for k in keys(kwargs) if k in atomic_symbols]

    A_X = format_A_X(metals_H, alpha_H, Dict(k => kwargs[k] for k in abunds_set))
    atm = interpolate_marcs(Teff, logg, A_X)

    synthesize_kwargs = [k => kwargs[k] for k in keys(kwargs) if !(k in abunds_set)]
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
    spectrum.wavelengths, flux, spectrum.cntm
end