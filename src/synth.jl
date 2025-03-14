"""
    synth(kwargs...)

This function creates a synthetic spectrum. It's easier to use than `synthesize`, but it gives you
less control. Unlike [`synthesize`](@ref), it **returns a tuple of `(wavelengths, rectified_flux, cntm)`**
(Wavelength in Å, rectified flux as a unitless number between 0 and 1, and continuum in erg/s/cm^5).

# Keyword arguments

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
               vmic=1.0,
               synthesize_kwargs=Dict(),
               format_A_X_kwargs=Dict(),
               abundances...,)
    A_X = format_A_X(metals_H, alpha_H, abundances; format_A_X_kwargs...)
    atm = interpolate_marcs(Teff, logg, A_X)

    # synthesize kwargs currently must be symbols, which is annoying
    spectrum = synthesize(atm, linelist, A_X, wavelengths; vmic, synthesize_kwargs...)
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