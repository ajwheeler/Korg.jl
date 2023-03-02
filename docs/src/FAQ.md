### What's the difference between `Korg.species` and `Korg.Species`?
The short answer is that they do the same thing, but the lower-case version can only be applied to string literals.
Both of these are ways to constructing objects of type [`Korg.Species`](@ref), which is used to represent atom and molecules with specific charges, e.g. C II or neutral FeH.  
`Korg.Species` is the constructor, which acts like a normal function.  `Korg.species` is a [non-standard string literal](https://docs.julialang.org/en/v1/manual/metaprogramming/#meta-non-standard-string-literals), which is used without parentheses, like this: `Korg.species"Mg I"`.
It can only be applied to string literals (things that you actually type out between quotes, not variables containing strings.)
The string macro does the work of constructing the species object at compile time, which can make code much faster in specific circumstances.
