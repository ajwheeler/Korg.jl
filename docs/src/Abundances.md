# Abundances

Korg's functions take arguments specifying abundances in two different ways.
[`Korg.synthesize`](@ref) takes a vector of ``A(X)`` abundances, one for every element from H to U.
([`Korg.format_A_X`](@ref) is a convienient way to create this vector.) Other functions
([`synth`](@ref), [`Korg.Fit.fit_spectrum`](@ref), etc.) take arguments specifying abundances more
indirectly. In either case, it's important to be clear about about how the arguments are
interpreted.

We'll focus on [`Korg.format_A_X`](@ref), which is what other functions use under the hood.
But first, some notation.
