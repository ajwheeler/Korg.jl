# [Wavelengths](@id wldocs)

Many of Korg's core functions operate on uniformly spaced wavelength sequences (note: all wavelengths have units of Å).
For every such such function the wavelength sequence(s) can be specified in any of the following ways:
- a pair of values `(λstart, λend)` specifies the bounds of wavelength sequence and a triple `(λstart, λend, λstep)` allows customization over spacing. `(λstart, λend)` is equivalent to `(λstart, λend, 0.01)`
- a vector of `N` pairs and/or triples, `(λstart1, λstop1[, λstep1]), ..., (λstartN, λstopN[, λstepN])` specifies `N` disjoint wavelength sequences.

## Examples

It is instructive to consider a few illustrative examples.
Although these examples use [`Korg.synth`](@ref), but apply equally to all functions that take wavelength
parameters, e.g. [`Korg.synthesize`](@ref), [`Korg.Fit.fit_spectrum`](@ref), [`Korg.RV_prec_from_noise`](@ref), and [`Korg.apply_LSF`](@ref).

First let's consider a simple continuous wavelength range:

```@example 1
using Korg, PythonPlot # hide
figure() # hide #md
wls, flux, _ = synth(Teff=5777, logg=4.44,
                    wavelengths=(5050, 5100)) # wavelength range (in Å)
plot(wls, flux)
ylabel("flux")
xlabel(L"$\lambda$ [\AA]")
gcf() # hide #md
```

In many cases, you want to compute spectra in noncontinuous chunks.
We can specify these chunks by passing a vector of pairs.

```@example 1
figure() # hide #md
wls, flux, _ = synth(Teff=5777, logg=4.44,
                     wavelengths=[(5050, 5080), (5090, 5100)])
plot(wls, flux)
ylabel("flux")
xlabel(L"$\lambda$ [\AA]")
gcf() # hide #md
```

The spacing for any wavelength sequence specified by a pair is 0.01 Å.
When you specify a sequence with a triple, the third element is used to overwrite the spacing with an arbitrary value.
Let's sample every 5 Å in the first window and every 3 Å in the second (while this sampling is too coarse for accurate spectra, these choices make it easy to visualize differences from the previous example).

```@example 1
figure() # hide #md
wls, flux, _ = synth(Teff=5777, logg=4.44,
                     wavelengths=[(5050, 5080, 5), (5110, 5120, 2)])
scatter(wls, flux)
ylabel("flux")
xlabel(L"$\lambda$ [\AA]")
gcf() # hide #md
```

## Gotchas and edge cases

- For a given tuple, `(λstart, λstop, λstep)`, `λend` may not always be included.
  For example, if `λstop` - `λstart` is not an integer multiple of `λstep`, the last last wavelength
  produced will be less than `λstop`.
- Wavelength sub-ranges _must not overlap_.  For example, `[(5000, 5010, 0.01), (5010, 5020, 0.2)]`,
  will fail, because `5010` is duplicated. Likewise, `[(5000, 5010, 0.02), (5009.9, 5020, 0.2)]`, will
  also fail
- Wavelenth sub-range must be sorted.  For example, `[(6000, 7000), (3000, 4000)]` will fail.
- All tuples must satisfy the invariants `1 < λstart < λstop` and `0 < λstep <= λstop - λstart`.
  Violations of these invariants will produce undefined behavior.
  - Historically, Korg inferred **(in an inconsistent manner)** that wavelengths have units of cm when `λstart<1.0`. You should **not** expect/rely upon this behavior.
  - The Korg developers reserve the right to freely modify the behavior when the invariants are violated, at any time (i.e. in a minor release).
- Because Julia's range objects are not guaranteed to produce bitwise identical ranges across Julia
  versions, Korg's wavelengths also lack this guarantee.

## Internals: [`Korg.Wavelengths`](@ref)

Functions that take wavelength parameters as arguments all pass them immediately to the
[`Korg.Wavelengths`](@ref) constructor, which turns them into the internal representation, so
the documentation for that function is also a good place to look if you are digging into the
Korg internals.
