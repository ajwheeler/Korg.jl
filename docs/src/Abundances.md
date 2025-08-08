# [Abundances](@id abundances)

Korg's functions take arguments specifying abundances in two different ways.
[`Korg.synthesize`](@ref) takes a vector of ``A(X)`` abundances, one for every element from H to U.
([`Korg.format_A_X`](@ref) is a convenient way to create this vector.) Other functions
([`synth`](@ref), [`Korg.Fit.fit_spectrum`](@ref), etc.) take arguments specifying abundances more
indirectly. In either case, it's important to be clear about how the arguments are
interpreted.

The "abundance" of an element, ``X``, can refer to both ``A(X)``
("absolute") and \[``X``/H\] ("solar relative"). These are defined as:

```math
A(X) = \log_{10} \left( \frac{n_X / n_\mathrm{H}} \right) + 12
```

and

```math
[X / \mathrm{H}] = A(X) - \left(A(X)\right)_\odot = \log_{10} \left ( \frac{n_X n^\odot_\mathrm{H}}{n_\mathrm{H} n^\odot_X} \right )
```

More fraught than "abundance" is "metallicity", which is routinely used to refer to the metal mass
fraction (``Z``), the iron abundance (solar relative or absolute), or to

```math
[M / \mathrm{H}]  = \log_{10} \left( \frac{\sum_\mathrm{X \in M} n_X}{n_\mathrm{H}} \frac{n^\odot_\mathrm{H}}{\sum_\mathrm{X \in M} n^\odot_X} \right),
```

where the set ``M`` can refer to everything heavier than He, or to a subset of elements that
researchers were able to measure. [`Korg.get_metals_H`](@ref) computes this.

```@example 1
using Korg # hide
Korg.get_metals_H(Korg.default_solar_abundances)
```

To generate abundances with a different metallicity, do:

```@example 1
A_X = Korg.format_A_X(-1)
Korg.get_metals_H(A_X)
```

But watch out! If we also specify \[N/H\], we get a different result:

```@example 1
A_X = Korg.format_A_X(-1, Dict(["N" => -0.5]))
Korg.get_metals_H(A_X)
```

That's because the first argument of [`Korg.format_A_X`](@ref) specifies what the "default"
abundance for each element should be, not the total metallcity.
It would be easy for it to scale the abundances of everything except N so that the resulting
metallicity was still `-1.0`, but this would be confusing in its own way.
