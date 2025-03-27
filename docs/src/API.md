This page documents the Korg API for users of the package. Low-level functions that only people
digging into the internals of Korg will be interested in can be found in the
[Developer documentation](@ref).

# [Top-level functions](@id API)
If you are synthesize a spectrum Korg, these are the functions you will call.
These functions are exported, so if you do `using Korg`, you can call them unqualified (i.e.
`synthesize` instead of `Korg.synthesize`).

```@docs
synth
synthesize
read_linelist
load_ExoMol_linelist
read_model_atmosphere
interpolate_marcs
format_A_X
```

## Built-in linelists
```@docs
Korg.get_APOGEE_DR17_linelist
Korg.get_GALAH_DR3_linelist
Korg.get_GES_linelist
Korg.get_VALD_solar_linelist
```

# Fitting
```@docs
Korg.Fit.fit_spectrum
Korg.Fit.ews_to_abundances
Korg.Fit.ews_to_stellar_parameters
```

# Secondary functions
You don't need use these to synthesize spectra, but they might be relevant depending on what you are
doing.

## Postprocessing
These are used to transform observational or synthetic spectra.

```@docs
Korg.apply_LSF
Korg.compute_LSF_matrix
Korg.air_to_vacuum
Korg.vacuum_to_air
```

## Misc

```@docs
Korg.blackbody
Korg.prune_linelist
Korg.merge_close_lines
```
