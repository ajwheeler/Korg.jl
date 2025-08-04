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
plot(λs, flux, "k-")
xlabel("λ [Å]")
ylabel("Flux")
gcf() # hide #md

# The [`synth`](@ref) function is a one-stop shop for synthesizing a spectrum.  It takes many keyword
# arguments, but the only required ones are `Teff` (the effective temperature in Kelvin) and `logg`
# (the log surface gravity in cgs units).  Korg will automatically construct the model atmosphere
# by interpolating from its built-in grid of MARCS model atmospheres (taking abundances, as well as
# `Teff` and `logg`, into account).
#
# ## Choose a linelist
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

vald_lines = Korg.get_VALD_solar_linelist() # lines for the Sun from 3000 to 9000 Å

# ## Specify the abundances
#
# To specify the abundances for a synthesis, the `M_H` and `alpha_H` keyword arguments specify the
# default metallicity and alpha enhancement, respectively.  To specify the hydrogen-relative,
# abundances of individual elements, pass their atomic symbol as a keyword argument to `synth`.
# There are unavoidable subtleties in abundance notation, so see section TODO for the messy details.
# For now, just know that the "metallicity" of a given mixture is not necessarily the same as the
# input `M_H` keyword argument, depending on how "metallicity" is defined.
#
# Let's synthesize another spectrum with the linelist we've selected, and with specific abundances.
# We'll also demonstrate a few more keyword arguments to [`synth`](@ref), but
# [`read the docs`](@ref synth) for more.

λs, flux, cntm = synth(; Teff=5777, logg=4.44, wavelengths=(5000, 5050), # what we had before
                       linelist=vald_lines, # use the linelist we selected
                       M_H=-1.1, alpha_H=-1.0, # metal-poor, mildy alpha-enhanced star
                       C=-0.5, # [C/H] = -0.5 dex (carbon enhanced)
                       R=20000, # Simulate and LSF with resolving power, R = λ/Δλ of 20,000
                       vsini=7, # projected rotational velocity of 7 km/s
                       vmic=2) # microturbulence of 2 km/s
plot(λs, flux, "k-")
xlabel("λ [Å]")
ylabel("Flux")
gcf() # hide #md

# # Going deeper: using `synthesize`
#
# The [`synth`](@ref) function provides a convenient interface that aims to cover most use cases,
# but in some circumstances you might want more control. The [`synthesize`](@ref) function is the
# more low-level function upon which [`synth`](@ref) is built. [`synthesize`](@ref) has a few
# required arguments.  Its signature is:
#
#     synthesize(atm, linelist, A_X, wavelenths...; kwargs...)
#
# [`synthesize`](@ref) takes a model atmosphere, a linelist, a vector of abundances, and parameter
# specifying the wavelengths to use as required arguments. We've covered linelists above, so let's
# take the rest in turn.
#
# ## Abundances
#
# For the sake of unambiguousness, [`synthesize`](@ref) takes a vector of abundances, `A_X`, as an
# for all elements from H to U. Note that unlike the [`synth`](@ref) function (by default), these
# are _not_ solar relative. They are in "standard" ``A(X)`` format,
# ``A(X) = \log_{10} \left( \frac{n_X}{n_H} \right)``.
# You can of course construct this vector yourself, but it's often easier to use the
# [`format_A_X`](@ref) function to do it for you.
# The interface for [`format_A_X`](@ref) is similar to that of [`synth`](@ref), but a tiny bit more
# explicit.  Instead of specifying abundances of individual elements with keyword arguments,
# we pass a dictionary of element symbols and abundances.

metal_poor_A_X = format_A_X(-0.5) # [M/H] = -1/2
alpha_rich_A_X = format_A_X(0, 0.5) # all [M/H] = 0, but [alpha/H] = 0.5]
Ni_enriched_A_X = format_A_X(Dict("Ni" => 1.0)) # all [M/H] = 0, except [Ni/H] = 1.0
display(Ni_enriched_A_X)
