# API

This page docuements the Korg API for users of the package. Low-level functions that only people 
digging into the internals of Korg will be interested in can be found on the 
[Developer documentation](@ref).

# Top-level functions
If you are synthesize a spectrum Korg, these are the functions you will call.  
These functions are exported, so if you do `using Korg`, you can call them unquallified (i.e.
`synthesize` instead of `Korg.synthesize`).  

```@docs
synthesize
read_linelist
read_model_atmosphere
interpolate_marcs
format_A_X
```

# Secondary functions
You don't need use these to synthesize spectra, but they might be relevant depending on what you are 
doing.

## Postprocessing
These are used to transform observational or synthetic spectra.

```@docs
Korg.rectify
Korg.constant_R_LSF
```

## Line absorption 
These functions can be used to directly compute line opacities. 

```@docs
Korg.line_absorption!
Korg.line_profile
Korg.hydrogen_line_absorption!
Korg.setup_hydrogen_stark_profiles
Korg.voigt
```

## Continuum absorption
These function can be used to directly compute continuum absorption coefficients.

##### Continuum Absorption Kwargs

All bound-free and free-free absorption functions (other than `absorption_coef_bf_TOPBase`) support
the following keyword arguments:

- `error_oobounds::Bool`: specifies the function's behavior when it encounters invalid `ν` or `T`
  values. When `false` (the default), the function assumes that the absorption coefficient for
  these invalid `ν` or `T` values is `0.0`. Otherwise, a `DomainError` is thrown when invalid `ν`
  or `T` values are encountered.
- `out_α::Union{Nothing,AbstractVector}`: When this is `nothing` (the default case), the function
  will simply allocate a new vector to store the output continuum absorption coefficients.
  Alternatively, this can be a vector (with the same length as the function's `ν` argument). In
  this case, the function will directly add the computed continuum absorption to the elements of
  this vector, in-place (the vector is also returned by the function).

```@docs
Korg.total_continuum_absorption
Korg.ContinuumAbsorption.H_I_bf
Korg.ContinuumAbsorption.H_I_ff
Korg.ContinuumAbsorption.Hminus_bf
Korg.ContinuumAbsorption.Hminus_ff
Korg.ContinuumAbsorption.H2plus_bf_and_ff
Korg.ContinuumAbsorption.He_II_bf
Korg.ContinuumAbsorption.He_II_ff
Korg.ContinuumAbsorption.Heminus_ff
Korg.ContinuumAbsorption.electron_scattering
Korg.ContinuumAbsorption.rayleigh
Korg.ContinuumAbsorption.hydrogenic_bf_absorption
Korg.ContinuumAbsorption.hydrogenic_ff_absorption
```

## Statistical mechanics
These functions can be used to calculate the number densities of all species in a given atospheric 
layer, or other statmech calculations. 

```@docs
Korg.saha_ion_weights
Korg.molecular_equilibrium_equations
Korg.molecular_equilibrium
```

## Misc

```@docs
Korg.blackbody
Korg.get_absolute_abundances
Korg.air_to_vacuum
Korg.vacuum_to_air
```
