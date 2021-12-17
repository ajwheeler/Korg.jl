This page contains info for people interested in contributing code to Korg.  If you 
have questions, do not hesitate to ask.

## Code guidelines
- Try to be explicit about units throughout the code, particularly when not using CGS.
- Whenever possible, calculations should be precise up to a factor of $$10^{-3}$$.  When it's easy and inexpensive, they should be precise to $$10^{-5}$$ or better.  
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

## Experimental support for bound-free metal absorption coefficients

### Overview

Korg currently has experimental support for computing the absorption coefficient contribution
from bound-free absorption by metals using data from the Opacity Project's online database,
[TOPbase](http://cdsweb.u-strasbg.fr/topbase/topbase.html).

For a given species (e.g. neutral C or singly-ionized C) TOPbase provides many tables of bound-free
cross-sections (without stimulated emission correction), as functions of wavelength. The 
`weighted_bf_cross_section_TOPBase` can be used to compute the effective bound-free cross-section
(including corrections from stimulated emission), assuming LTE, for the given species.

Our tentative plan for the future is to use `weighted_bf_cross_section_TOPBase` for each species
to construct a 2D tables of the effective bound-free cross-section that can be interpolated over
temperature and wavelength. However, in the short-term the `absorption_coef_bf_TOPBase` function
(which calls `weighted_bf_cross_section_TOPBase`) can be used to directly compute the absorption
coefficient from the TOPbase data.

### Data Sources

Currently, to call `weighted_bf_cross_section_TOPBase` or `absorption_coef_bf_TOPBase` for a
given species, a table is required from TOPbase that is the result of a
[photoionisation Cross Section query](http://cdsweb.u-strasbg.fr/topbase/xsections.html).
The query should only request data for a single species.

Precise details about the requirements of the query is provided in the docstring of
`each_photo_ion_subtable`. Note that leading/trailing white-space can cause issues for the
parsers of these tables. For examples of these tables, see the `TOPbase_cross_section_He_II.txt`
and `TOPbase_cross_section_H_I.txt` files in the `test/data` directory.

The path to this table should be directly passed to `weighted_bf_cross_section_TOPBase` and
`absorption_coef_bf_TOPBase` via the `cross_sec_file` keyword argument. Alternatively, the code
supports searching for the table at the path `$KORG_OP_DIR/<symbol>_<ion-state>.txt`, where
- `KORG_OP_DIR` is the name of an environment variable
- `<symbol>` is replaced with the species atomic symbol
- `<ion-state>` is replaced with the capitalized Roman numerals describing the ionization state.
Examples of the base filename (under this environment variable approach) include `H_I.txt`,
`He_I.txt`, `He_II.txt`, or `Al_III.txt`.

In the future, we also plan to add support for reading in photoionization tables from the
[NORAD-Atomic-Data webpage](https://norad.astronomy.osu.edu/#codes). 

!!! warning "Caveats"

    Be aware that there are discrepancies between the theoretical energy levels computed by the
    Opacity Project and empirical measurements. This can be seen by looking at the plots produced
    by the `misc/gen_Opacity_Project_comparison_plot.jl` plot. To recover accurate UV opacities,
    [Gustafsson+ 2008](https://ui.adsabs.harvard.edu/abs/2008A%26A...486..951G/abstract) needed to
    correct the precise photon energies in the tables for C I, Mg I, Al I, and Si I. Before
    providing full support for these opacities, we may want to consider making similar corrections.


## Complete API
Here are all the documented methods in Korg.

```@autodocs
Modules = [Korg, Korg.ContinuumAbsorption]
```
