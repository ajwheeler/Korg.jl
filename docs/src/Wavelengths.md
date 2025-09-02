# [Wavelengths](@id wldocs)

Many of Korg's core functions operate on uniformly spaced wavelength sequences (note: all wavelengths have units of Å).
For every such such function expects the wavelength sequence(s) can be specified in any of the following ways:
- a pair of values `(λstart, λend)` specifies the bounds of wavelength sequence and a triple `(λstart, λend, λstep)` allows customization over spacing. `(λstart, λend)` is equivalent to `(λstart, λend, 0.01)`
- a vector of `N` pairs and/or triples, `(λstart1, λstop1[, λstep1]), ..., (λstartN, λstopN[, λstepN])` specifies `N` disjoint wavelength sequences.

## Examples

It is instructive to consider a few illustrative examples.
Although these examples use [`Korg.synth`](@ref), but apply equally to all functions that take wavelength
parameters, e.g. [`Korg.synthesize`](@ref), [`Korg.Fit.fit_spectrum`](@ref), [`Korg.RV_prec_from_noise`](@ref), and [`Korg.apply_LSF`](@ref).

First let's consider a simple continuous wavelength range:

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

The spacing for any wavelength sequence specified by a pair is 0.01 Å.
When you specify a sequence with a triple, the third element is used to overwrite the spacing with an arbitrary value.
Let's sample every 5 Å in the first window and every 3 Å in the second (while this sampling is too coarse for accurate spectra, these choices make it easy to visualize differences from the previous example).

```@example 1
wls, flux, _ = synth(Teff=5777, logg=4.44,
                     wavelengths=[(5050, 5080, 5), (5090, 5100, 3)])
scatter(wls, flux)
ylabel("flux")
xlabel(L"$\lambda$ [\AA]")
gcf() # hide
```


## Internals: [`Korg.Wavelengths`](@ref)

Functions that take wavelength parameters as arguments all pass them immediately to the
[`Korg.Wavelengths`](@ref) constructor, which turns them into the internal representation, so
the documentation for that function is also a good place to look if you are digging into the
Korg internals.
