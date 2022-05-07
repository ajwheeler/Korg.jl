using Interpolations: LinearInterpolation, Throw
using HDF5: h5read

using ..ContinuumAbsorption: hydrogenic_bf_absorption, hydrogenic_ff_absorption, ionization_energies

const _H_I_ion_energy = ionization_energies[1][1] # not sure if this is a good idea
const _H⁻_ion_energy = 0.7552 # eV

_H_I_bf(ν, T, nH_I_div_partition, ion_energy = _H_I_ion_energy, nmax_explicit_sum = 8,
        integrate_high_n = true) =
           hydrogenic_bf_absorption(ν, T, 1, nH_I_div_partition, ion_energy, nmax_explicit_sum,
                                    integrate_high_n)
"""
    H_I_bf(ν, T, nH_I_div_partition, [ion_energy], [nmax_explicit_sum], [integrate_high_n];
           kwargs...)

Compute the bound-free linear absorption coefficient contributed by all energy states of a
neutral Hydrogen atom.

# Required Arguments
- `ν::AbstractVector{<:Real}`: sorted frequency vector in Hz
- `T`: temperature in K
- `nH_I_div_partition` is the total number density of neutral Hydrogen divided by the its
   partition function.

# Optional Arguments
- `ion_energy`: The ionization energy of Hydrogen. By default, this is set to the values loaded 
  into the global `ionization_energies` list.
- `nmax_explicit_sum`: The highest energy level whose absorption contribution is included
   in the explicit sum. The contributions from higher levels are included in an integral.
- `integrate_high_n::bool`: When this is `false`, bf absorption from higher energy states are not
   estimated at all. Default is `true`.

For a description of the kwargs, see [Continuum Absorption Kwargs](@ref).

# Notes
This function wraps [`hydrogenic_ff_absorption`](@ref). See that function for implementation
details.
"""
H_I_bf = bounds_checked_absorption(_H_I_bf; ν_bound = Interval(0, Inf),
                                   temp_bound = Interval(0, Inf))

_H_I_ff(ν, T, nH_II, ne) = hydrogenic_ff_absorption(ν, T, 1, nH_II, ne)

"""
    H_I_ff(ν, T, nH_II, ne; kwargs...)

Compute the H I free-free linear absorption coefficient α.

The naming scheme for free-free absorption is counter-inutitive. This actually refers to the
reaction:  `photon + e⁻ + H II -> e⁻ + H II`.

# Arguments
- `ν::AbstractVector{<:Real}`: sorted frequency vector in Hz
- `T`: temperature in K
- `nH_II`: the number density of ionized Hydrogen in cm⁻³.
- `ne`: the number density of free electrons.

For a description of the kwargs, see [Continuum Absorption Kwargs](@ref).

# Notes
This function wraps [`hydrogenic_ff_absorption`](@ref). See that function for implementation
details.
"""
H_I_ff = bounds_checked_absorption(_H_I_ff; ν_bound = λ_to_ν_bound(_gauntff_λ_bounds),
                                   temp_bound = _gauntff_T_bounds)

"""
    _ndens_Hminus(nH_I_div_partition, ne, T, ion_energy = _H⁻_ion_energy)

Compute the number density of H⁻ (implements eqn 5.10 of Kurucz 1970). This is an application of
the saha equation where the "ground state" is H⁻ and the "first ionization state" is H I. The
partition function of H⁻ is 1 at all temperatures.
"""
function _ndens_Hminus(nH_I_div_partition, ne, T, ion_energy = _H⁻_ion_energy)
    # The value of _H_I_ion_energy energy doesn't matter at all becasue the exponential evaluates
    # to 1 for n = 1
    nHI_groundstate = ndens_state_hydrogenic(1, nH_I_div_partition, T, _H_I_ion_energy)

    # should should probably factor this out
    if T < 1000
        # this isn't a hard requirement at all. This is just a sanity check in case an argument got
        # forgotten (which is allowed since ion_energy has a default value)
        # In reality, we should use the lower temperature limit in which the exponential evaluates
        # to Inf. That should probably be our domain limit
        throw(DomainError(T, "Unexpected T value."))
    end

    # coef = (h^2/(2*π*m))^1.5
    coef = 3.31283018e-22 # cm³*eV^1.5
    β = 1.0/(kboltz_eV*T)
    0.25 * nHI_groundstate * ne * coef * β^1.5 * exp(ion_energy * β)
end

const _Hminus_bf_cross_section = let
    fn = joinpath(_data_dir, "McLaughlin2017Hminusbf.h5")
    ν = h5read(fn, "nu")
    σ = h5read(fn, "sigma")
    LinearInterpolation(ν, σ, extrapolation_bc=Throw())
end

function _Hminus_bf(ν::Real, T::Real, nH_I_div_partition::Real, ne::Real)
    cross_section = _Hminus_bf_cross_section(ν) # in units of cm²
    stimulated_emission_correction = (1 - exp(-hplanck_cgs*ν/(kboltz_cgs*T)))
    n_Hminus = _ndens_Hminus(nH_I_div_partition, ne, T, _H⁻_ion_energy)
    n_Hminus * cross_section * stimulated_emission_correction
end

"""
    Hminus_bf(ν, T, nH_I_div_partition, ne; kwargs...)

Compute the H⁻ bound-free linear absorption coefficient α

# Arguments
- `ν::AbstractVector{<:Real}`: sorted frequency vector in Hz
- `T`: temperature in K
- `nH_I_div_partition`: the total number density of H I divided by its partition function.
- `ne`: the electron number density

For a description of the kwargs, see [Continuum Absorption Kwargs](@ref).

# Notes
This function assumes that n(H⁻) ≪ n(H I) + n(H II). The number density of n(H⁻) isn't precomputed 
as part of Korg's molecular equlibrium, it's computed here instead.

This function is adapted from equation 8.12 from Grey (2005). That equation specifies the
absorption coefficient per H I atom (uncorrected by stimulated emission) as:

``8.316 \\times 10^{-10} \\sigma_{bf}(H^-) P_e \\theta^{2.5} 10^{(\\theta \\chi_{\\mathrm{H}^-})} U(H⁻,T)/ U(H I,T)``

where:
- ``\\sigma_{bf}(\\mathrm{H}^-)`` is the photo dissociation cross-section (Grey uses the variable
  ``\\alpha_{bf}`` in place of ``\\sigma_{bf}``). We estimate this by linearly interpolating data
  from [McLaughlin (2017)](https://ui.adsabs.harvard.edu/abs/2017JPhB...50k4001M/abstract).
- ``\\chi_{\\mathrm{H}^-}`` is the ionization enrgy of H⁻.
- θ = log10(e)/(k*T) or θ = 5040/T in units of eV⁻¹
- U(H⁻,T) is the partition function of H⁻ at temperature T. This is always 1
- U(H I,T) is the partition function of H I at temperature T. This is 2 at low to intermediate T.
Eqn 8.12 of Grey (2005) implicitly assumes that (U(H⁻,T)/ U(H I,T)) is always equal to 0.5.

This expression is simply a rewritten form of ``\\sigma_{bf}(H^-) n(H^-)/n(H I),`` where
``n(H^-)/n(H I)`` has been replaced with the expanded form of the saha equation, in which
n₀ = n(H⁻) and n₁ = n(H I).

Equations 8.18 and 8.19 of Gray (2005), indicate that the H⁻ bound-free opacity (with the
stimulated emission correction) is given by:

``\\kappa_\\nu = \\sigma_{bf}(H^-) \\frac{n(H^-)}{n(H I)} (1 - \\exp \\left( \\frac{-h\\nu}{k T}\\right)) \\frac{n(H I)}{n(H I) + n(H II)} \\frac{n(H I) + n(H II)}{\\rho}``

In other words, the linear absorption coefficient is: ``\\alpha_\\nu = \\sigma_{bf}(H^-) n(H⁻) (1 - \\exp \\left( \\frac{-h\\nu}{k T}\\right))``

Wishart (1979) expects the tabulated data to have better than 1% percent accuracy. Mathisen (1984)
suggests that this data has better than 3% accuracy.
"""
Hminus_bf = bounds_checked_absorption(
    _Hminus_bf;
    ν_bound = closed_interval(1.8238892857120884e14, 2.417989242625068e19),
    temp_bound = Interval(0, Inf)
)

function _Hminus_ff(ν::Real, T::Real, nH_I_div_partition::Real, ne::Real)
    λ = c_cgs*1e8/ν # in Angstroms

    logλ = log10(λ)
    log2λ = logλ * logλ
    log3λ = log2λ * logλ
    log4λ = log3λ * logλ

    f0 =  -2.2763 -   1.6850 * logλ +  0.76661 * log2λ -  0.053346 * log3λ
    f1 = +15.2827 -   9.2846 * logλ +  1.99381 * log2λ -  0.142631 * log3λ
    f2 = -197.789 + 190.266  * logλ - 67.9775  * log2λ + 10.6913   * log3λ - 0.625151 * log4λ

    logθ = log10(5040.0/T)
    αff_H⁻ = 1e-26 * 10.0^(f0 + f1 * logθ + f2 * (logθ * logθ))
    # Pₑ * α_ff(H⁻) gives the absorption coefficient in units of cm² per ground state H I atom
    Pₑ = ne * kboltz_cgs * T

    # Account for the fact that n(H I, n=1) might be slightly smaller than the entire number
    # density of H I. There is only really a difference at the highest temperatures. For the
    # temperature range where this approximation is valid, less than 0.23% of all H I atoms are not
    # in the ground state.

    # this calculation could reduce to nHI_gs = 2.0*nH_I_div_partition
    nHI_gs = ndens_state_hydrogenic(1, nH_I_div_partition, T, _H_I_ion_energy)

    return αff_H⁻ * Pₑ * nHI_gs
end

"""
    Hminus_ff(ν, T, nH_I_div_partition, ne; kwargs...)

Compute the H⁻ free-free linear absorption coefficient α

The naming scheme for free-free absorption is counter-inutitive. This actually refers to the
reaction:  photon + e⁻ + H I -> e⁻ + H I.

# Arguments
- `ν::AbstractVector{<:Real}`: sorted frequency vector in Hz
- `T`: temperature in K
- `nH_I_div_partition::Flt`: the total number density of H I divided by its partition function.
- `ne`: the number density of free electrons.

For a description of the kwargs, see [Continuum Absorption Kwargs](@ref).

# Notes
This is taken from equation 8.13 of Gray (2005). The equation uses a polynomial fig against Table 1
of [Bell & Berrington's (1987)](https://doi.org/10.1088/0022-3700/20/4/019), which tabulates values
for "the H⁻ absorption coefficient" (the values include the correction for stimulated emission). We
denote this quantity as `f(H⁻)` (Grey instead denotes it with the variable, ``\\alpha_{ff}(H^-)``).
The polynomial fits `f(H⁻)` in the range 2520 K ≤ T ≤ 10080 K and 2604 Å ≤ λ ≤ 113918 Å.
According to Grey, the polynomial fit to the tabulated data typically has 1% precision. We find
that at worst, the discrepancy never exceeds 2.25%.

Based on equations 8.13, 8.18, and 8.19 of Gray (2005), the free-free linear absorption coefficient
of H⁻ is given by: `α = f(H⁻) * Pₑ * n(H I)`.

Based on Section 5.3 from Kurucz (1970), I'm fairly confident that the "more correct" version of
this equation should actually read `α = f(H⁻) * Pₑ * n(H I, n = 1)`, and that the form inferred
from Gray (2005) implicitly assumes that `n(H I, n = 1) ≈ n(H I)`. In detail, the Boltzmann
equation relates `n(H I, n = 1)` and `n(H I)` as:
```
 n(H I, n = 1) = n(H I)*gₙ₌₁/U(T)*exp(-Eₙ₌₁/(k*T)) = n(H I) * 2/U(T)*exp(0) = n(H I) * 2/U(T).
 ```
Based on the values of the partition function in Table D.2 of Grey, this approximation seems fairly
good. It only introduces a ≤0.5% overestimate in n(H I, n = 1) over the temperature range where the
polynomial is valid.

We also considered the polynomial fit in Section 5.3 from Kurucz (1970). Unfortunately, it seems
to be wrong (it gives lots of negative numbers).
"""
Hminus_ff = bounds_checked_absorption(_Hminus_ff,
                                      ν_bound = λ_to_ν_bound(closed_interval(2.604e-5,1.13918e-3)),
                                      temp_bound = closed_interval(2520, 10080))

function _H2plus_bf_and_ff(ν::Real, T::Real, nH_I::Real, nH_II::Real)
    λ = c_cgs*1e8/ν # in ångstroms
    β_eV = 1.0/(kboltz_eV * T)
    stimulated_emission_correction = (1 - exp(-hplanck_eV*ν*β_eV))

    K = Stancil1994.K_H2plus(T) # n(H I) n(H II) / n(H₂⁺)
    σbf = Stancil1994.σ_H2plus_bf(λ, T)
    σff = Stancil1994.σ_H2plus_ff(λ, T)

    (σbf/K + σff) * nH_I * nH_II * stimulated_emission_correction
end

"""
    H2plus_bf_and_ff(ν, T, nH_I_div_partition, n_HII; kwargs...)

Compute the combined H₂⁺ bound-free and free-free linear absorption coefficient α using the tables 
from [Stancil 1994](https://ui.adsabs.harvard.edu/abs/1994ApJ...430..360S/abstract).

# Arguments
- `ν::AbstractVector{<:Real}`: sorted frequency vector in Hz
- `T`: temperature in K
- `nH_I`: the total number density of H I divided by its partition 
   function.
- `nH_II`: the number density of H II (not of H₂⁺).

For a description of the kwargs, see [Continuum Absorption Kwargs](@ref).

# Notes
This computes the H₂⁺ number density from those of H I and H II, since ionized molecules 
are not included in the molecular equlibrium calculation.  Stancil provides approximate equilibrium 
constants, K, for the molecule, but the Barklem and Collet values used elsewhere by Korg may be more 
reliable.  Once ionized molecules are fully supported, those values should be used instead. While 
cross sections are tabulated down to only 3150 K, the cross sections could be linearly interpolated 
1000 K or so lower, if reliable K values are availble.

Because n(H₂⁺) is computed on the fly, the "cross-sections" used (internally) by this function have units of 
cm^-5, since they must be multiplied by n(H I) and n(H II) to obtain absorption coefficients.

Stancil, [Bates (1952)](https://ui.adsabs.harvard.edu/abs/1952MNRAS.112...40B/abstract), and
Gray (2005) all use n(H I), but Kurucz (1970) notes in section 5.2 that the photodisociation of H₂⁺ 
produces a ground state H atom and that we should therefore use n(H I, n=1) instead. This difference 
causes a descrepancy of a fraction of a percent (at most) for T ≤ 1.2e4 K. Here we use n(H I).
"""
H2plus_bf_and_ff = bounds_checked_absorption(
    _H2plus_bf_and_ff;
    ν_bound = λ_to_ν_bound(closed_interval(7e-6, 2e-3)),
    temp_bound = closed_interval(3150, 25200)
)

