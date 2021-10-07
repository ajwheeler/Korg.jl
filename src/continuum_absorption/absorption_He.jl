"""
According to Gray (2005), the bound-free contributions from He⁻ are usually assumed to be 
negligible because it only has one bound level with an ionization energy 19 eV. Supposedly the 
population of that level is too small to be worth considering.

We are currently missing free-free and bound free contributions from He I.
"""

using ..ContinuumAbsorption: hydrogenic_bf_absorption, hydrogenic_ff_absorption, ionization_energies
const _He_II_ion_energy = ionization_energies[2][2] # not sure if this is a good idea

_He_II_bf(ν, T, nHe_II_div_partition, ion_energy = _He_II_ion_energy, nmax_explicit_sum = 9,
          integrate_high_n = true) =
              hydrogenic_bf_absorption(ν, T, 2, nHe_II_div_partition, ion_energy, nmax_explicit_sum,
                                       integrate_high_n)
"""
    He_II_bf(ν, T, nHe_I_div_partition, [ion_energy], [nmax_explicit_sum], [integrate_high_n])

Compute the bound-free linear absorption coefficient contributed by all energy states of a
singly-ionized Helium atom.

# Required Arguments
- `ν`: frequency in Hz
- `T`: temperature in K
- `nHe_II_div_partition` is the total number density of singly-ionized Helium divided by the its
   partition function.

# Optional Arguments
- `ion_energy`: The ionization energy of Helium. By default, this is set to the values loaded 
  into the global `ionization_energies` list.
- `nmax_explicit_sum`: The highest energy level whose absorption contribution is included
   in the explicit sum. The contributions from higher levels are included in an integral.
- `integrate_high_n::bool`: When this is `false`, bf absorption from higher energy states are not
   estimated at all. Default is `true`.

# Notes
This function simply wraps [`hydrogenic_ff_absorption`](@ref). See that function for additional
details.
"""
He_II_bf = bounds_checked_absorption(_He_II_bf; ν_bound = Interval(0, Inf),
                                     temp_bound = Interval(0, Inf))

_He_II_ff(ν, T, nHe_III, ne) = hydrogenic_ff_absorption(ν, T, 2, nHe_III, ne)

"""
    He_II_ff(ν, T, nHe_III, ne)

Compute the He II free-free linear absorption coefficient α.

The naming scheme for free-free absorption is counter-inutitive. This actually refers to the
reaction:  `photon + e⁻ + He III -> e⁻ + He III`.

#Arguments
- `ν`: frequency in Hz
- `T`: temperature in K
- `nH_III`: the number density of doubly-ionized Helium in cm⁻³.
- `ne`: the number density of free electrons.


# Notes
This function wraps [`hydrogenic_ff_absorption`](@ref). See that function for implementation
details.
"""
He_II_ff = bounds_checked_absorption(
    _He_II_ff,
    ν_div_T_bound = closed_interval( (1e-4, 10^1.5)./(hplanck_eV/kboltz_eV)... ),
    temp_bound = closed_interval( (4*RydbergH_eV/kboltz_eV)./(1e2, 1e-3)... )
)

# Compute the number density of atoms in different He I states
# taken from section 5.5 of Kurucz (1970)
function ndens_state_He_I(n::Integer, nsdens_div_partition::Real, T::Real)
    g_n, energy_level = if n == 1
        (1.0, 0.0)
    elseif n == 2
        (3.0, 19.819)
    elseif n == 3
        (1.0, 20.615)
    elseif n == 4
        (9.0, 20.964)
    else
        # we definitely could add more energy levels
        throw(DomainError(n," Unknown excited state properties for He I"))
    end
    nsdens_div_partition * g_n * exp(-energy_level/(kboltz_eV *T))
end

function _Heminus_ff(ν::Real, T::Real, nHe_I_div_partition::Real, ne::Real)

    λ = c_cgs * 1.0e8 / ν # Å
    if !(5063.0 <= λ <= 151878.0)
        throw(DomainError(λ, "The wavelength must lie in the interval [5063 Å, 151878 Å]"))
    elseif !(2520.0 <= T <= 10080.0)
        throw(DomainError(T, "The temperature must lie in the interval [2520 K, 10080 K]"))
    end

    θ = 5040.0 / T
    θ2 = θ * θ
    θ3 = θ2 * θ
    θ4 = θ3 * θ

    c0 =   9.66736 - 71.76242*θ + 105.29576*θ2 - 56.49259*θ3 + 10.69206*θ4
    c1 = -10.50614 + 48.28802*θ -  70.43363*θ2 + 37.80099*θ3 -  7.15445*θ4
    c2 =   2.74020 - 10.62144*θ +  15.50518*θ2 -  8.33845*θ3 +  1.57960*θ4
    c3 =  -0.19923 +  0.77485*θ -   1.13200*θ2 +  0.60994*θ3 -  0.11564*θ4

    logλ = log10(λ)

    # f includes contribution from stimulated emission
    f = 1e-26*10.0^(c0 + c1 * logλ + c2 * (logλ * logλ) + c3 * (logλ * logλ * logλ))

    Pe = ne*kboltz_cgs*T # partial pressure contributed by electrons
    nHe_I_gs = ndens_state_He_I(1, nHe_I_div_partition, T)

    f * nHe_I_gs * Pe
end

"""
    Heminus_ff(ν, T, nHe_I_div_partition, ne)

Compute the He⁻ free-free opacity κ.

The naming scheme for free-free absorption is counter-inutitive. This actually refers to the
reaction:  `photon + e⁻ + He I -> e⁻ + He I.`

# Arguments
- `ν`: frequency in Hz
- `T`: temperature in K
- `nHe_I_div_partition`: the total number density of H I divided by its partition function.
- `ne`: the number density of free electrons.

# Notes

This follows equation 8.16 of Grey (2005) which provides a polynomial fit absorption to data
tabulated in [John (1994)](https://ui.adsabs.harvard.edu/abs/1994MNRAS.269..871J/abstract).
According to that equation, the "continuous absorption coefficient per Hydrogen" is given by:
```
    f(He⁻_ff)*Pₑ*A(He) / (1 + Φ(He)/Pₑ)
```
in which
- `log_10(f(He⁻_ff))` is the polynomial term (Grey actually uses the variable α in place of f)
- A(He) is the number abundance of Helium particles relative to Hydrogen particles. For reference,
  ```
  A(He) = [n(He I) + n(He II) + n(He III)] / [n(H I) + n(H II)]
  ```
- Pₑ is the partial pressure contributed by electrons. For reference, `Pₑ = nₑ*kb * T`.
- `1/(1 + Φ(He)/Pₑ)` comes from the Saha equation and expresses `n(He I) / [n(He I) + n(He II)]`.

In the above expression, f(He⁻_ff)*Pₑ specifies the free-free atomic absorption coefficient per
ground state He I atom. The expression seems to implicitly assume that
- approximately all He I is in the ground state
- `n(He I) / [n(He I) + n(He II)]` is roughly `n(He I) / [n(He I) + n(He II) + n(HeIII)]`
(at least for the range of temperatures where f can be accurately computed).

To convert the expression to opacity, Grey divided it by the product of the mean molecular weight
and the Hydrogen mass. With that in mind, we can write the equation for linear absorption
coefficient, α_ν, (removing the apparent assumptions) as `α_ν = f(He⁻_ff)*Pₑ*n(He I, n=1).`

For 5063 Å ≤ λ ≤ 151878 Å and 2520 K ≤ T ≤ 10080 K, Gray (2005) claims that the polynomial fit the
tabulated data at a precision of 1% or better. In practice, we found that it only fits the data to
better than 3.1% (it's possible that for smaller λ the fit may be better). For reference, the
tabulated data in these ranges of values consist of an irregularly spaced rectangular grid with 15
λ values and 9 Temperature values. According to John (1994), improved calculations are unlikely to
alter the tabulated data for λ > 1e4Å, "by more than about 2%." The errors introduced by the
approximations for 5.06e3 Å ≤ λ ≤ 1e4 Å "are expected to be well below 10%."

An alternative approach using a fit to older data is provided in section 5.7 of Kurucz (1970).
"""
Heminus_ff = bounds_checked_absorption(
    _Heminus_ff,
    ν_bound = _λ_to_ν_bound(closed_interval(5.063e-5, 1.518780e-03)),
    temp_bound = closed_interval(2520, 10080)
)
