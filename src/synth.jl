"""
    synth(kwargs...)

This function creates a synthetic spectrum. It's easier to use than `synthesize`, but it gives you
less control. Rather than returning a `SynthesisResult`, it returns a tuple of
`(wavelengths, flux, cntm)`.

# Keyword arguments

TODO

# Returns

TODO
"""
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
               synthesize_kwargs=Dict(),
               abundances...,)
    # TODO prefer "X_H" to "X" for abundances
    A_X = format_A_X(metals_H, alpha_H, abundances)
    atm = interpolate_marcs(Teff, logg, A_X)

    # synthesize kwargs currently must be symbols, which is annoying
    spectrum = synthesize(atm, linelist, A_X, wavelengths; synthesize_kwargs...)
    flux = if rectify
        spectrum.flux ./ spectrum.cntm
    else
        spectrum.flux
    end
    if isfinite(R)
        # TODO this induces edge effects
        flux = apply_LSF(flux, spectrum.wavelengths, R)
    end
    if vsini > 0
        flux = apply_rotation(flux, spectrum.wavelengths, vsini)
    end

    spectrum.wavelengths, flux, spectrum.cntm
end