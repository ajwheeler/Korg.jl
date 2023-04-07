This pages lists major changes and additions to Korg since the publication of [the code paper](https://arxiv.org/abs/2211.00029). For a more detailed changelog, see [the github releases](https://github.com/ajwheeler/Korg.jl/releases).

- **0.17.0** introduced the use of the Mihalas-Daeppen-Hummer formalism to account for plasma effects when calculating Hydrogen lines and
  H I bound-free absorption.

- **0.15.0** added 28 new polyatomic molecules from exomol to the equilibrium calculations.

- **0.14.0** introduced a built-in model atmosphere interpolator, based on the SDSS MARCS grid.  See [`interpolate_marcs`](@ref).

- **0.10.0** is the version of Korg described by the code paper.
