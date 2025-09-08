This page documents the public, stable Korg API intended for users of the package.
Starting with the version 1.0 release, the developers commit to maintaining forward compatability for all of these functions (until, at least, the next major release).

If you are interested in contributing to Korg, we provide a more complete list of functions (that includes private functions) in the [Developer documentation](@ref).

# [Top-level functions](@id API)
If you are trying to synthesize a spectrum with Korg, these are the functions you will call.
These functions are exported, so if you do `using Korg`, you can call them unqualified (i.e.
`synthesize` instead of `Korg.synthesize`).

```@docs
synth
synthesize
format_A_X
interpolate_marcs
Korg.read_model_atmosphere
```

## Linelists
```@docs
Korg.read_linelist
Korg.load_ExoMol_linelist
Korg.get_APOGEE_DR17_linelist
Korg.get_GALAH_DR3_linelist
Korg.get_GES_linelist
Korg.get_VALD_solar_linelist
Korg.save_linelist
Korg.Line
```

# Fitting
```@docs
Korg.Fit.fit_spectrum
Korg.Fit.ews_to_abundances
Korg.Fit.ews_to_stellar_parameters
Korg.Fit.ews_to_stellar_parameters_direct
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
