# Developer documentation
This page contains info for people interested in contributing code to Korg.  If you 
have questions, do not hesitate to ask.

## Julia development
If you've never developed a package in Julia, here are some tips.
- Once you have a local copy of Korg, you can make it importable locally by `dev`ing it from the Julia environment in which you would like to run it: `Pkg.dev("/path/to/Korg")`.  If you run `Pkg.dev("Korg")` instead, Julia will automatically clone the repo to `~/.julia/dev/Korg`.  (See the [Pkg documentation](https://pkgdocs.julialang.org/) for details.)
- When working on your local copy of Korg, [Revise](https://github.com/timholy/Revise.jl) is easiest way to make and test changes without constantly restarting your Julia session.
- To run the test suite locally, start a Julia session in the Korg root directory and run `]activate .` then `test`.  This will automatically run `test/runtests.jl` in the test environment.
- Documentation is generated with [Documentor](https://github.com/JuliaDocs/Documenter.jl). Continuous integration on github will ensure that the documentation is generated without errors, but it won't catch all formatting problems.  If you wish to generate documentation locally in order to check that everything is as expected, start a Julia session in `Korg/docs`, activate the test environment (`]activate .`), and run `instantiate`.  This downloads and installs the docs dependences, and will only have to be run once.  To generate documentation, run `julia --project make.jl` on the command line (the `--project` flag activates the local environment).  The generated docs can be served from `docs/build`.

## Code guidelines
- Try to be explicit about units throughout the code, particularly when not using CGS.
- Whenever possible, calculations should be precise up to a factor of $$10^{-3}$$.  When it's easy and inexpensive, they should be precise to $$10^{-5}$$ or better.  There are exceptions to this (e.g. the Voigt function), but ideally they will eventually be vanquished.
- Ensure types are generic enough to support dual numbers and autodifferentiation. 
- Limit lines to 100 characters.
- Use ASCII characters in the names of functions that are part of the public API.
- Unless they will never be called elsewhere, provide docstrings describing the inputs, assumptions and outputs of any functions you write.

## Continuum absorption
Steps for implementing new continuum sources of absorption:
- Define a helper function that computes a single absorption coefficient (in units of cm⁻²). The function should accept `ν` (in Hz) and `T` (in K) as the first and second arguments, respectively. The convention is for it should share a name with the corresponding public function, but have an underscore pre-appended (e.g. we define `_H_I_bf` to help implement `H_I_bf`).
- The public function is the function constructed and returned by `Korg.ContinuumAbsorption.bounds_checked_absorption` that wraps the above helper function. This wrapper function implements bounds-checking for `ν` and `T` and supports the keyword arguments described in [Continuum Absorption Kwargs](@ref).
- Add a docstring describing the new function. At the very least, please describe any non-standard arguments and include a reference in the docstring to the source where the function was taken from.
- Add a line to `doc/src/API.md` under the `Continuum absorption` heading to render the docstring of your new function.
- Add a line to `total_continuum_absorption` that calls the new public function for absorption.

The first two steps may not apply for sources that don't directly depend on `ν` and `T` (e.g.
absorption from scattering).

## `Manifest.toml` and `Project.toml`
[Read more about these files here.](https://pkgdocs.julialang.org/v1/toml-files/)
These files are used by the julia package manager.  `Project.toml` records dependencies, and you'll 
notice that the test and docs directories have their own `Project.toml`s for test and 
documentation-specific dependencies.  `Manifest.toml` records exact package versions used when a 
package was run.  It enables someone else to reproduce the exact environment and results later.  
There are a few directories containing scripts that generate Korg's data files which have there own 
`Projects.toml`s and `Manifest.toml`s, for example `data/bf_cross-sections/`.

## Where to put data 
If you are adding data to Korg, `data/README` provides an overview of the options and how to decide 
between them.

## Complete API
Here are all the documented methods in Korg.

```@autodocs
Modules = [Korg, Korg.Fit, Korg.CubicSplines, Korg.ContinuumAbsorption, Korg.ContinuumAbsorption.Stancil1994,
Korg.ContinuumAbsorption.Peach1970, Korg.RadiativeTransfer, Korg.RadiativeTransfer.MoogStyleTransfer, Korg.RadiativeTransfer.BezierTransfer]
```
