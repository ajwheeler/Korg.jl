# Wavelengths

Many Korg functions takes wavelengths as one of their parameters.  Because Korg allows for
flexible wavelength specification, it's helpful to go over the possibilities.


First is the simple case: a continuous wavelength range:

```@example 1
using Korg, PythonPlot # hide
wls, flux, _ = synth(Teff=5777, logg=4.44,
                    wavelengths=(5050, 5100)) # wavelength range (in Å)
plot(wls, flux)
ylabel("flux")
xlabel(L"$\lambda$ [\AA]")
gcf() # hide
```

In many cases, you want to compute spectra in noncontinuous chunks.  To do that, instead of passing
a pair, pass a vector of pairs.

```@example 1
wls, flux, _ = synth(Teff=5777, logg=4.44,
                     wavelengths=[(5050, 5080), (5090, 5100)])
scatter(wls, flux)
ylabel("flux")
xlabel(L"$\lambda$ [\AA]")
gcf() # hide
```

By default, Korg uses samples wavelength every 0.01 Å, but you can configure this by specifying a
different value as the last wavelength parameter. Let's sample every 5 Å in the first window and
every 3 Å in the second (this is too coarse for accurate spectra, but it makes it easy to see what's
going on).

```@example 1
wls, flux, _ = synth(Teff=5777, logg=4.44,
                     wavelengths=[(5050, 5080, 5), (5090, 5100, 3)])
scatter(wls, flux)
ylabel("flux")
xlabel(L"$\lambda$ [\AA]")
gcf() # hide
```


!!! note
      Functions that take wavelength parameters as arguments all pass them immediately to the
      [`Korg.Wavelengths`](@ref) constructor, which turns them into the internal representation, so
      the documentation for that function is also a good place to look.

These examples use [`Korg.synth`](@ref), but apply equally to all functions that take wavelength
parameters, e.g. [`Korg.synthesize`](@ref), [`Korg.Fit.fit_spectrum`](@ref),
[`Korg.RV_prec_from_noise`](@ref), and [`Korg.apply_LSF`](@ref).
