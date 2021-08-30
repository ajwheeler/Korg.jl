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
``weighted_bf_cross_section_TOPBase`` can be used to compute the effective bound-free cross-section 
(including corrections from stimulated emission), assuming LTE, for the given species.

Our tentative plan for the future is to use ``weighted_bf_cross_section_TOPBase`` for each species
to construct a 2D tables of the effective bound-free cross-section that can be interpolated over
temperature and wavelength. However, in the short-term the ``absorption_coef_bf_TOPBase`` function
can be used to directly compute the absorption coefficient from the TOPbase data.

### Data Sources

Currently, to call ``weighted_bf_cross_section_TOPBase`` or ``absorption_coef_bf_TOPBase`` for a
given species, 2 tables are required from TOPbase that are the result of
1. a [photoionisation Cross Section query](http://cdsweb.u-strasbg.fr/topbase/xsections.html) of
   the given species.
2. an [energy levels query](http://cdsweb.u-strasbg.fr/topbase/energy.html) that includes the
   target species (this query can include data for multiple species of a given element)

Precise details about the requirements of each query are respectively provided in the docstrings of
`each_photo_ion_subtable` and `read_electron_configurations`. Note that leading/trailing
white-space can cause issues for the parsers of these tables. For examples of these tables, see the
``TOPbase_cross_section_He_II.txt`` and ``TOPbase_electron_config_He.txt`` files in the
``test/data`` directory.

The paths to these tables should be directly passed to ``weighted_bf_cross_section_TOPBase`` and
``absorption_coef_bf_TOPBase`` via the ``cross_sec_file`` and ``elec_conf_file`` keyword arguments.
Alternatively, the ``KORG_OP_CROSS_SECTION_DIR`` and ``KORG_OP_ELECTRON_CONFIG_DIR`` environment
variables can be used to specify different directories that hold these files. When using the
environment variable approach:
- cross-section files must be located at ``$KORG_OP_CROSS_SECTION_DIR/<symbol>_<ion-state>.txt``,
  where ``symbol`` is replaced with the species atomic symbol and <ion-state> is replaced with
  the upper-case capital letters describing the ionization state. Examples of the base filename
  include ``H_I.txt``, ``He_I.txt``, or ``He_II.txt``
- electron-configuration files must be located at ``$KORG_OP_ELECTRON_CONFIG_DIR/<symbol>.txt``.
  Note, that in this case properties for all relevant species of a given element should be
  specified in a single file. For example, if you want to separately compute absorption for
  ``He_I`` and ``He_II``, the states for both species (when using the environment variables) are
  expected to be saved to a single file called ``He.txt``.

We plan to drop the need for the electron-configuration tables in the future as it provides
redundant information (it's currently only used for legacy reasons). In the future, we also plan to
add support for reading in photoionization tables from the
[NORAD-Atomic-Data webpage](https://norad.astronomy.osu.edu/#codes).

### Caveats

Be aware that there are discrepancies between the theoretical energy levels computed by the Opacity
Project and empirical measurements. This can be seen by looking at the plots produced by the
``misc/gen_Opacity_Project_comparison_plot.jl`` plot. To recover accurate UV opacities,
[Gustafsson+ 2008](https://ui.adsabs.harvard.edu/abs/2008A%26A...486..951G/abstract) needed to
correct the precise photon energies in the tables for C I, Mg I, Al I, and Si I. Before providing
full support for these opacities, we may want to consider making similar corrections.


## Complete API
Here are all the documented methods in Korg.

```@autodocs
Modules = [Korg, Korg.ContinuumOpacity]
```
