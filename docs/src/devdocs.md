This page contains info for people interested in contributing code to Korg.  If you 
have questions, do not hesitate to ask.

## code guidelines
- Try to be explicit about units throughout the code, particularly when not using CGS.
- Whenever possible, calculations should be precise up to a factor of $$10^{-3}$$.  When it's easy and inexpensive, they should be precise to $$10^{-5}$$ or better.  
- Ensure types are generic enough to support dual numbers and autodifferentiation. 
- Limit lines to 100 characters.

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
Modules = [Korg, Korg.ContinuumOpacity]
```
