## What's the difference between `Korg.species` and `Korg.Species`?
Short answer: you can always use `Korg.Species` like so: `Korg.Species("C II")`.
Longer answer: they do the same thing, but the lower-case version can only be applied to string literals.
Both of these are ways to constructing objects of type [`Korg.Species`](@ref), which is used to represent atom and molecules with specific charges, e.g. C II or neutral FeH.
`Korg.Species` is the constructor, which acts like a normal function.  `Korg.species` is a [non-standard string literal](https://docs.julialang.org/en/v1/manual/metaprogramming/#meta-non-standard-string-literals), which is used without parentheses, like this: `Korg.species"Mg I"`.
It can only be applied to string literals (things that you actually type out between quotes, not variables containing strings.)
The string macro does the work of constructing the species object at compile time, which can make code much faster in specific circumstances.

## Are Korg spectra bitwise reproducible?
Korg spectra should be bitwise reproducible under these conditions:
- Using a single thread (Julia's default configuration)
- Running on the same machine
- Using identical Korg and dependency versions
- Providing identical input parameters

## What's the deal with Korg's flux units?
[`Korg.synthesize`](@ref) and [`Korg.synth`](@ref) (when set to not rectify) return fluxes in units of `erg/s/cm^-4/Å`. The "extra" factor of `cm^-2` is included to account for the projected area of the star, which Korg makes no assumptions about.  To get the usual astrophysical flux, multiply by ``R_\mathrm{star}^2``.

When running with multiple threads, small numerical differences may occur due to floating-point operations being performed in different orders.
