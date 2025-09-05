# [Download this page as a Jupyter notebook](./Basics.ipynb) #md
# You may have to execute these cells to see the plot outputs #nb

using Suppressor # hide #md
@suppress begin # hide #md
    using Korg, PythonPlot
end # hide #md

# When you import with `using`, a few functions are imported directly into the namespace
# ([`synth`](@ref), [`synthesize`](@ref), [`interpolate_marcs`](@ref), [`format_A_X`](@ref)).  All
# other functions need to be accessed with the `Korg` prefix, e.g. [`Korg.air_to_vacuum`](@ref).
# These examples use the PythonPlot package, which provides a nice Julia interface to Matplotlib,
# but of course you can plot however you like.

# # Synthesizing a spectrum the easy way: `synth`
# Synthesizing a spectrum is easy!

λs, flux, cntm = synth(; Teff=5777, logg=4.44, wavelengths=(5000, 5050))
figure(; figsize=(12, 4)) # hide #md
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
# Julia is just-in-time compiled, which means that the first time in each session you call a method
# (e.g. `synth`), the code will be compiled, which can take several seconds.  The first time you call a
# method, it may be slow, with lots of time devoted to compilation. When you call it a second time,
# it's much faster, because no code needs to be compiled. Note that no data is being cached and
# re-used here. Try calling `synth` again (changing the parameters if you like) to see the speedup.
#
# You probably want to be able to control more than just the effective temperature, surface gravity,
# and wavelength range. Let's go over some of the other things you can set when calling
# [`synth`](@ref).
#
# ## Choose a linelist
#
# A "linelist" is a list of atomic and molecular transitions to be included in the synthetic
# spectrum.  Korg supports the VALD, MOOG, Kurucz, and ExoMol format linelists (see
# [`Korg.read_linelist`](@ref), [`Korg.load_ExoMol_linelist`](@ref)).
# It also has several built-in linelists for convenience:
#    - [`Korg.get_VALD_solar_linelist`](@ref) for the Sun from 3000 to 9000 Å.
#    - [`Korg.get_APOGEE_DR17_linelist`](@ref) for the latest APOGEE linelist
#    - [`Korg.get_GALAH_DR3_linelist`](@ref) for the GALAH DR3 and DR4 linelist
#    - [`Korg.get_GES_linelist`](@ref) for the Gaia ESO survey linelist
#
# Each linelist is a vector of `Korg.Line` objects, but you don't need to worry about the details if you
# just want to pass it into Korg.
#

vald_lines = Korg.get_VALD_solar_linelist() # lines for the Sun from 3000 to 9000 Å

# ## Specifying abundances
#
# To specify the abundances for a synthesis, the `M_H` and `alpha_H` keyword arguments specify the
# default metallicity and alpha enhancement, respectively.  To specify the hydrogen-relative
# abundances of individual elements, pass their atomic symbol as a keyword argument to `synth`.
# There are unavoidable subtleties in abundance notation, so see [Abundances](@ref) for the messy
# details.
# For now, just know that the "metallicity" of a given mixture is not necessarily the same as the
# input `M_H` keyword argument, depending on how "metallicity" is defined.
#
# ##  A more complex example #nb
# ## [A more complex example](@id synth-example) #md
#
# Let's synthesize another spectrum with the linelist we've selected, and with specific abundances.
# We'll also demonstrate a few more keyword arguments to [`synth`](@ref), but
# you can [read the docs](@ref synth) for more.

λs, flux, cntm = synth(; Teff=5777, logg=4.44, wavelengths=(5000, 5050), # what we had before
                       linelist=vald_lines, # use the linelist we selected
                       M_H=-1.1, alpha_H=-1.0, # metal-poor, mildly alpha-enhanced star
                       C=-0.5, # [C/H] = -0.5 dex (carbon enhanced)
                       R=20000, # Simulate an LSF with resolving power, R = λ/Δλ of 20,000
                       vsini=7, # projected rotational velocity of 7 km/s
                       vmic=2) # microturbulence of 2 km/s
figure(; figsize=(12, 4)) # hide #md
plot(λs, flux, "k-")
xlabel("λ [Å]")
ylabel("Flux")
gcf() # hide #md

# # Going deeper: `synthesize`
#
# The [`synth`](@ref) function provides a convenient interface that aims to cover most use cases,
# but in some circumstances you might want more control. The [`synthesize`](@ref) function is the
# more low-level function upon which [`synth`](@ref) is built. [`synthesize`](@ref) has a few
# required arguments.  Its signature is:
#
#     synthesize(atm, linelist, A_X, wavelengths...; kwargs...)
#
# [`synthesize`](@ref) takes a model atmosphere, a linelist, a vector of abundances, and parameters
# specifying the wavelengths to use as required arguments. We've covered linelists above, so let's
# take abundances and model atmospheres in turn.
#
# ## Abundances
#
# For the sake of unambiguity, [`synthesize`](@ref) takes a vector of abundances, `A_X`, as an
# argument for all elements from H to U. Note that unlike the [`synth`](@ref) function (by default), these
# are _not_ solar relative. They are in "standard" ``A(X)`` format,
# ``A(X) = \log_{10} \left( \frac{n_X}{n_H} \right)``.
# You can of course construct this vector yourself, but it's often easier to use the
# [`format_A_X`](@ref) function to do it for you.
# The interface for [`format_A_X`](@ref) is similar to that of [`synth`](@ref), but a tiny bit more
# explicit.  Instead of specifying abundances of individual elements with keyword arguments,
# we pass a dictionary of element symbols and abundances.

metal_poor_A_X = format_A_X(-0.5) # [M/H] = -1/2
alpha_rich_A_X = format_A_X(0, 0.5) # all [M/H] = 0, but [alpha/H] = 0.5
Ni_enriched_A_X = format_A_X(Dict("Ni" => 1.0)) # all [M/H] = 0, except [Ni/H] = 1.0

# ## Model Atmosphere
#
# As mentioned above, Korg can interpolate from a grid of MARCS model atmospheres that it is shipped
# with using the [`interpolate_marcs`](@ref) function. These are primarily
# [those generated for the Sloan Digital Sky Survey (SDSS)](https://www.sdss4.org/dr17/irspec/apogee-libraries/),
# but the grid has been augmented with more atmospheres at lower metallicities.
# See [the docs](@ref interpolate_marcs) for more details on the models and how they get
# interpolated. We can create a model atmosphere object to pass to [`synthesize`](@ref) like this:

atm = interpolate_marcs(5777, 4.44, Ni_enriched_A_X) # solar Teff and logg, Ni-enriched abundances

# Korg will look at your abundance vector and translate it into the appropriate parameters for the
# atmosphere grid.  It's possible to bypass this and specify the abundance parameters directly, but
# it's easy to mess up and not recommended.
#
# If you have a model atmosphere file you want to use instead of Korg's internal grid, you can read
# it in with [`Korg.read_model_atmosphere`](@ref).

# ## Putting it all together
#
# Now we can synthesize a stellar spectrum by passing the linelist and atmosphere to
# [`synthesize`](@ref), along with upper and lower wavelengths (in Å). Korg uses wavelengths *in
# vacuo*, but you can use [`Korg.air_to_vacuum`](@ref) and [`Korg.vacuum_to_air`](@ref) to convert
# back and forth.

res = synthesize(atm, vald_lines, format_A_X(0), (4000, 4015));

# The object returned by [`synthesize`](@ref) is a [`Korg.SynthesisResult`](@ref), which contains
# lots of information. As with [`synth`](@ref), the wavelengths, flux, and continuum are available,
# though note that in a
# [`Korg.SynthesisResult`](@ref), the flux has not had the continuum divided out (i.e. the spectrum
# has not been continuum normalized, or rectified).

figure(; figsize=(12, 4)) # hide #md
plot(res.wavelengths, res.flux, "k-") # "./" is how you do element-wise division in Julia
xlabel(L"$\lambda$ [Å]")
ylabel("flux [erg/s/cm^2/Å]")
gcf() # hide #md

# ## Per-species number densities
#
# One of the pieces of information in a [`Korg.SynthesisResult`](@ref) is the per-species number
# densities, which are stored as a dictionary mapping [`Korg.Species`](@ref) objects to arrays
# with one entry for each layer in the model atmosphere.

figure() # hide #md
temps = Korg.get_temps(atm)
# These strings represent different species. Below, we pass them to Korg.Species to construct a Species object.
for spec in ["H I", "H II", "O I", "OH"]
    plot(temps, res.number_densities[Korg.Species(spec)]; label=spec)
end

legend()
yscale("log")
xlabel(L"$T$ [K]")
ylabel(L"$n$ [cm$^{-3}$]")
gcf() # hide #md

# ## The absorption coefficient
# Another piece of information in a [`Korg.SynthesisResult`](@ref) is the absorption coefficient,
# α, in units of cm^-1.

figure(; figsize=(12, 4)) # hide #md
plot(res.wavelengths, res.alpha'); # res.alpha' is the adjoint of res.alpha.
yscale("log")
xlabel(L"$\lambda$ [Å]")
ylabel(L"$\alpha$ [cm$^{-1}$]");
gcf() # hide #md

# # Post-processing your spectrum
#
# Given a synthesized spectrum, you may want to manipulate it further to make it look more
# like real data. Korg provides functions for applying a line-spread function (LSF) and rotational
# broadening to a spectrum. Note that [the example above](@ref synth-example) demonstrates shortcuts
# to do this with [`synth`](@ref).
#
# ## [Line-spread function (LSF)](@id LSF)
# Korg provides two ways of applying an LSF to a spectrum. The first, [`Korg.apply_LSF`](@ref),
# is simpler (this is what [`synth`](@ref) does). It returns a new flux vector and does not modify
# the wavelength sampling.

R = 10_000 # the resolving power, R = λ/Δλ
low_res_flux = Korg.apply_LSF(res.flux, res.wavelengths, R)

figure(; figsize=(12, 4)) # hide #md
plot(res.wavelengths, res.flux, "k-"; label="high-res")
plot(res.wavelengths, low_res_flux, "C1-"; label="low-res")
legend()
xlabel(L"$\lambda$ [Å]")
ylabel("flux [erg/s/cm^2/Å]")
gcf() # hide #md

# The second way of applying an LSF is slightly more powerful: [`Korg.compute_LSF_matrix`](@ref).
# This function computes an efficiently-represented sparse matrix that transforms an
# infinite-resolution flux vector to a low-resolution flux vector, simultaneously resampling the
# wavelengths. This saves a lot of computing time if you need to apply the same LSF to many spectra.

new_wavelengths = 4000:0.1:4015
LSF = Korg.compute_LSF_matrix(res.wavelengths, new_wavelengths, R)

# The `LSF` object is a sparse matrix that we can use to transform the high-resolution flux vector
# (or any other computed on the same wavelength grid)

resampled_low_res_flux = LSF * res.flux

figure(; figsize=(12, 4)) # hide #md
plot(res.wavelengths, low_res_flux, "C1-"; label="low-res")
plot(new_wavelengths, resampled_low_res_flux, "C2-"; label="low-res (resampled wavelengths)")
legend()
xlabel(L"$\lambda$ [Å]")
ylabel("flux [erg/s/cm^4/Å]")
gcf() # hide #md

# ### Varying LSFs
#
# If your LSF has a Gaussian width that varies with wavelength, that's no problem!  Just pass a
# function that maps from wavelength (in Å) to R.

resolution(λ) = 10_000 + 2000 * sin(λ) # R varying with wavelength
crazy_lsf_flux = Korg.apply_LSF(res.flux, res.wavelengths, resolution)

figure(; figsize=(12, 4)) # hide #md
plot(res.wavelengths, res.flux, "k-"; label="high-res")
plot(res.wavelengths, low_res_flux, "C1-"; label="low-res (constant LSF)")
plot(res.wavelengths, crazy_lsf_flux, "C3-"; label="low-res (varying LSF)")
legend()
xlabel(L"$\lambda$ [Å]")
ylabel("flux [erg/s/cm^4/Å]")
gcf() # hide #md

# ## Rotational broadening
#
# To apply the broadening resulting from differential line-of-sight velocity across the face of the
# star due to rotation, you can use [`Korg.apply_rotation`](@ref).

vsini = 30 # km/s
rot_flux = Korg.apply_rotation(res.flux, res.wavelengths, 15)

figure(; figsize=(12, 4)) # hide #md
plot(new_wavelengths, LSF * res.flux, "C1-"; label="original")
plot(new_wavelengths, LSF * rot_flux, "C4-"; label="with rotation")
legend()
xlabel(L"$\lambda$ [Å]")
ylabel("flux [erg/s/cm^4/Å]")
gcf() # hide #md
