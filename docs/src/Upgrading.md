
# Upgrading to Korg `v1.0`

This page explains the changes in Korg v1.0 that might trip you up if you are coming from an earlier version.

- [`synthesize`](@ref), [`Korg.MolecularCrossSection`](@ref), and [`Korg.prune_linelist`](@ref) now always
  take only a single argument for wavelength specification.
  Bounds can be passed as tuple, e.g. `synthesize(atm, linelist, A_X, (λ_lower, λ_upper))`.
  For `synthesize`, passing more than one argument will trigger a deprecation warning for now.
- Keyword arguments named `m_H` are now called `M_H`.  This is more common notation, and is meant to
  look less like it means "hydrogen mass". This affects [`synth`](@ref),
  [`Korg.Fit.ews_to_stellar_parameters`](@ref), and [`Korg.Fit.fit_spectrum`](@ref)
- The default solar abundances ([`Korg.default_solar_abundances`](@ref)) are now
  _Bergemann et al. 2025_ ([`Korg.bergemann_2025_solar_abundances`](@ref)).
- The output units of flux are now `erg/s/cm^4/Å`, not `erg/s/cm^5`, which is more consistent with
  Korg's wavelength units.
- [`synth`](@ref) now requires that you specify at least `Teff` and `logg`.  All other arguments are
  still optional.
- When synthesizing in the infrared, by default, the Mihalas-Daeppen-Hummer formalism is no longer
  used to self-consistently adjust hydrogen level populations for the purpose of bound-bound
  transitions. Previously, this was a warning.
- `air_wavelengths` is no longer an allowed keyword argument for [`Korg.synthesize`](@ref). Users are
  recommended to use a workflow that involves explicit conversions with [`Korg.air_to_vacuum`](@ref) and [`Korg.vacuum_to_air`](@ref).
  The old functionality is still available in [`Korg.Wavelengths`](@ref), but this workflow can be
  confusing and we recommend avoiding it.
- The (previously deprecated) `air_wavelengths` keyword argument for [`Korg.synthesize`](@ref) has been removed.
