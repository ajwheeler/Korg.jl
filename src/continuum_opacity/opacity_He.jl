"""
According to Gray (2005), the bound-free contributions from He⁻ are usually assumed to be 
negligible because it only has one bound level with an ionization energy 19 eV. Supposedly the 
population of that level is too small to be worth considering.

We are currently missing free-free and bound free contributions from He I.
"""

using ..ContinuumOpacity: hydrogenic_bf_opacity, hydrogenic_ff_opacity, ionization_energies
const _He_II_ion_energy = ionization_energies["He"][2] # not sure if this is a good idea

He_II_bf(nHe_II_div_partition, ν, ρ, T, ion_energy = _He_II_ion_energy,
         nmax_explicit_sum = 9, integrate_high_n = true) =
             hydrogenic_bf_opacity(1, nHe_II_div_partition, ν, ρ, T, ion_energy, nmax_explicit_sum,
                                   integrate_high_n)
# He II free-free actually refers to the reaction: photon + e⁻ + He III -> e⁻ + He III.
He_II_ff(nHe_III, ne, ν, ρ, T) = hydrogenic_ff_opacity(2, nHe_III, ne, ν, ρ, T)


# Compute the number density of atoms in different He I states
# taken from section 5.5 of Kurucz (1970)
function ndens_state_He_I(n::Integer, nsdens_div_partition::Flt, T::Flt) where {Flt<:AbstractFloat}
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


"""
    Heminus_ff(nHe_I_div_partition, ne, ν, ρ, T)

Compute the He⁻ free-free opacity κ.

The naming scheme for free-free absorption is counter-inutitive. This actually refers to the
reaction:  photon + e⁻ + He I -> e⁻ + He I.

# Arguments
- `nHe_I_div_partition::Flt`: the total number density of H I divided by its partition function.
- `ne::Flt`: the number density of free electrons.
- `ν::Flt`: frequency in Hz
- `ρ::Flt`: mass density in g/cm³
- `T::Flt`: temperature in K

# Notes

This follows equation 8.16 of Grey (2005) which provides a polynomial fit absorption to data
tabulated in [John 1994](https://ui.adsabs.harvard.edu/abs/1994MNRAS.269..871J/abstract). According
to that equation the "continuous absorption coefficient per Hydrogen" is given by:
    α(He⁻_ff)*Pₑ*A(He) / (1 + Φ(He)/Pₑ)
in which
- log_10(α(He⁻_ff)) is the polynomial term.
- A(He) is the number abundance of Helium particles relative to Hydrogen particles. For reference,
  A(He) = [n(He I) + n(He II) + n(He III)] / [n(H I) + n(H II)]
- Pₑ is the partial pressure contributed by electrons. For reference, Pₑ = nₑ*kb * T.
- 1/(1 + Φ(He)/Pₑ) comes from the Saha equation and expresses n(He I) / [n(He I) + n(He II)].

In the above expression, α(He⁻_ff)*Pₑ specifies the free-free atomic absorption coefficient per
ground state He I atom. The expression seems to implicitly assume that
- approximately all He I is in the ground state
- n(He I) / [n(He I) + n(He II)] is roughly n(He I) / [n(He I) + n(He II) + n(HeIII)]
(at least for the range of temperatures where α can be accurately computed).

To convert the expression to opacity, they divided it by the product of the mean molecular weight
and the Hydrogen mass. With that in mind, we can write an the equation for opacity (removing the
apparent assumptions) as κ_ν = α(He⁻_ff)*Pₑ*n(He I, n=1)/ρ.

For 5063 Å ≤ λ ≤ 151878 Å and 2520 K ≤ T ≤ 10080 K, Gray (2005) claims that the polynomial fit the
tabulated data at a precision of 1% or better. In practice, we found that it only fits the data to
better than 3.1% (it's possible that for smaller λ the fit may be better). For reference, the
tabulated data in these ranges of values consist of an irregularly spaced rectangular grid with 15
λ values and 9 Temperature values. According to John (1994), improved calculations are unlikely to
alter the tabulated data for λ > 1e4Å, "by more than about 2%." The errors introduced by the
approximations for 5.06e3 Å ≤ λ ≤ 1e4 Å "are expected to be well below 10%."

An alternative approach using a fit to older data is provided in section 5.7 of Kurucz (1970).
"""
function Heminus_ff(nHe_I_div_partition::Flt, ne::Flt, ν::Flt, ρ::Flt,
                    T::Flt) where {Flt<:AbstractFloat}

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

    # α includes contribution from stimulated emission
    α = 1e-26*10.0^(c0 + c1 * logλ + c2 * (logλ * logλ) + c3 * (logλ * logλ * logλ))

    Pe = ne*kboltz_cgs*T # partial pressure contributed by electrons
    nHe_I_gs = ndens_state_He_I(1, nHe_I_div_partition, T)

    α * nHe_I_gs * Pe / ρ
end
