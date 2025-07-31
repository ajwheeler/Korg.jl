# When you import with "using", a few functions are imported directly into the namespace (`synth`,
# `synthesize`, `interpolate_marcs`, `format_A_X`).  All other need to be accessed with the `Korg`
# prefix, e.g. `Korg.air_to_vacuum`.
# We are also using the PythonPlot package, which provides a nice Julia interface to Matplotlib, 
# but of course you can plot however you like.

using Suppressor # hide #md
@suppress begin # hide #md
    using Korg, PyPlot
end # hide #md

# # Synthesizing a spectrum the easy way: `synth`
# Synthesizing a spectrum is easy!

λs, flux, cntm = synth(; Teff=5777, logg=4.44, wavelengths=(5000, 5050))
plot(λs, flux)
xlabel("λ [Å]")
ylabel("Flux")
gcf() # hide #md

# # Choose a linelist
# 
# A "linelist" is a list of atomic and molecular transitions to be included in the synthetic
# spectrum.  Korg supports the VALD, MOOG, Kurucz, and ExoMol format linelists (see 
# [`Korg.read_linelist`](@ref), [`Korg.load_ExoMol_linelist`](@ref)). 
# It also has several built-in for convenience:
#    - [`Korg.get_VALD_solar_linelist`](@ref) for the Sun from 3000 to 9000 Å
#    - [`Korg.get_APOGEE_DR17_linelist`](@ref) for latest APOGEE linelist
#    - [`Korg.get_GALAH_DR3_linelist`](@ref) for the GALAH DR3 and DR4 linelist
#    - [`Korg.get_GES_linelist`](@ref) for the Gaia ESO survey linelist
#
# Each linelist is a vector of `Korg.Line` objects, but you don't need to worry about details if you
# just want to pass it into Korg.
#

lines = Korg.get_VALD_solar_linelist() # lines for the Sun from 3000 to 9000 Å
