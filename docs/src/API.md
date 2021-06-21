This page docuements the Korg API for users of the package. Low-level functions that only people 
digging into the guts of Korg will be interested in can be found on the 
[developer documentation](../devdocs).

## top-level functions
If you are synthesize a spectrum Korg, these are the functions you will call.  
These functions are exported, so if you do `using Korg`, you can call them unquallified (i.e.
`synthesize` instead of `Korg.synthesize`).  

```@docs
synthesize
constant_R_LSF
read_line_list
read_model_atmosphere
```

## line absorption 
These functions can be used to directly compute line opacities. 

```@docs
Korg.line_absorption
Korg.line_profile
Korg.voigt
```

## continuum absorption
These function can be used to directly compute continuum opacities.

```@docs
Korg.total_continuum_opacity
Korg.ContinuumOpacity.H_I_bf
Korg.ContinuumOpacity.H_I_ff
Korg.ContinuumOpacity.Hminus_bf
Korg.ContinuumOpacity.Hminus_ff
Korg.ContinuumOpacity.H2plus_bf_and_ff
Korg.ContinuumOpacity.He_II_bf
Korg.ContinuumOpacity.He_II_ff
Korg.ContinuumOpacity.Heminus_ff
Korg.ContinuumOpacity.electron_scattering
```

## statistical mechanics
These functions can be used to calculate the number densities of all species in a given atospheric 
layer, or other statmech calculations. 

```@docs
Korg.saha_ion_weights
Korg.molecular_equilibrium_equations
Korg.molecular_equilibrium
```

## misc

```@docs
Korg.blackbody
Korg.get_absolute_abundances
```

