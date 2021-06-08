# we are missing Rayleigh Scattering

using ..ContinuumOpacity: hydrogenic_bf_opacity, hydrogenic_ff_opacity, ionization_energies

const _H_I_ion_energy = ionization_energies["H"][1] # not sure if this is a good idea
const _H⁻_ion_energy = 0.7552 # eV


H_I_bf(nH_I_div_partition, ν, ρ, T, ion_energy = _H_I_ion_energy,
       nmax_explicit_sum = 8, integrate_high_n = true) =
           hydrogenic_bf_opacity(1, nH_I_div_partition, ν, ρ, T, ion_energy, nmax_explicit_sum,
                                 integrate_high_n)
# H I free-free actually refers to the reaction: photon + e⁻ + H II -> e⁻ + H II.
H_I_ff(nH_II, ne, ν, ρ, T) = hydrogenic_ff_opacity(1, nH_II, ne, ν, ρ, T)

# compute the number density of H⁻ (implements eqn 5.10 of Kurucz 1970). This formula comes from
# inverting the saha equation, where n(H⁻) is n₀ and n(H I) is n₁. Note that U₀ = 1 at all
# temperatures.
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


"""
    _Hminus_bf_cross_section(ν)

Compute the H⁻ bound-free cross-section, which has units of cm^2 per H⁻ particle.

The cross-section does not include a correction for stimulated emission.
"""
function _Hminus_bf_cross_section(ν::Real)
    λ = c_cgs*1e8/ν # in Angstroms
    # we need to somehow factor out this bounds checking
    if !(2250 <= λ <= 15000.0)
        throw(DomainError(λ, "The wavelength must lie in the interval [2250 Å, 15000 Å]"))
    end

    λ2 = λ*λ
    λ3 = λ*λ2
    λ4 = λ*λ3
    λ5 = λ*λ4
    λ6 = λ*λ5

    αbf_H⁻ = (1.99654 - 1.18267e-5 * λ + 2.64243e-6 * λ2 - 4.40524e-10 * λ3 + 3.23992e-14 * λ4
              - 1.39568e-18 * λ5 + 2.78701e-23 * λ6)
    αbf_H⁻ * 1e-18
end


"""
    Hminus_bf(nH_I_div_partition, ne, ν, ρ, T, [ion_energy_H⁻])

Compute the H⁻ bound-free opacity κ

# Arguments
- `nH_I_div_partition::Flt`: the total number density of H I divided by its partition function.
- `ne` the electron number density
- `ν::Flt`: frequency in Hz
- `ρ::Flt`: mass density in g/cm³
- `T::Flt`: temperature in K
- `ion_energy_H⁻::Flt`: Specifies the ionization energy of the single state of H⁻ in eV. This is
   roughly 0.7552 eV.

# Notes
This function assumes that n(H⁻) ≪ n(H I) + n(H II). The number density of n(H⁻) should not be
pre-computed (instead it's computed internally by this function).

This function is adapted from equation 8.12 from Grey (2005). This equation gives the absorption
coefficient per H I atom (uncorrected by stimulated emission) is given by:
    αff_H⁻ * 8.316e-10 * Pₑ * θ^2.5 * 10^(ion_energy_H⁻ * θ) * (U(H⁻,T)/ U(H I,T))
where:
- αff_H⁻ is the photo dissociation cross-section. This can estimated with equation 8.11.
- θ = log10(e)/(k*T) or θ = 5040/T in units of eV⁻¹
- U(H⁻,T) is the partition function of H⁻ at temperature T. This is always 1
- U(H I,T) is the partition function of H I at temperature T. This is 2 at low to intermediate T.
Eqn 8.12 of Grey (2005) implicitly assumes that (U(H⁻,T)/ U(H I,T)) is always equal to 0.5.

This expression is simplicitly a rewritten form of: αff_H⁻ * n(H⁻)/n(H I) where n(H⁻)/n(H I) has
been replaced with the expanded form of the saha equation, in which n₀ = n(H⁻) and n₁ = n(H I).

Combining 8.18, and 8.19 of Gray (2005), indicate that the version opacity contribution of H⁻
bound-free absorption (with stimulated emission correction) is given by:

                      n(H⁻)                                  n(H I)          n(H I) + n(H II)
    κ_ν = α_bf(H⁻) * ------  * (1 - exp(-h*ν/(k*T))) * ------------------ * -----------------
                     n(H I)                             n(H I) + n(H II)            ρ

This can be rewritten as: κ_ν = α_bf(H⁻) * n(H⁻) * (1 - exp(-h*ν/(k*T))) / ρ

This function uses the polynomial provided in equation 8.11 of Gray (2005), that fits the tabulated
data from Wishart (1979). While Gray (2005) claims that the polynomial fits the data with 0.2%
precision for 2250 Å ≤ λ ≤ 15000 Å, in practice we find that it fits the data to better than 0.25%
precision. Wishart (1979) expects the tabulated data to have better than 1% percent accuracy.
"""
function Hminus_bf(nH_I_div_partition::Real, ne::Real, ν::Real, ρ::Real, T::Real, 
                   ion_energy_H⁻::Real = 0.7552)
    αbf_H⁻ = _Hminus_bf_cross_section(ν) # does not include contributions from stimulated emission
    stimulated_emission_correction = (1 - exp(-hplanck_cgs*ν/(kboltz_cgs*T)))
    n_H⁻ = _ndens_Hminus(nH_I_div_partition, ne, T, _H⁻_ion_energy)
    αbf_H⁻ * n_H⁻ * stimulated_emission_correction / ρ
end


"""
    Hminus_ff(nH_I_div_partition, ne, ν, ρ, T)

Compute the H⁻ free-free opacity κ

The naming scheme for free-free absorption is counter-inutitive. This actually refers to the 
reaction:  photon + e⁻ + H I -> e⁻ + H I.

# Arguments
- `nH_I_div_partition::Flt`: the total number density of H I divided by its partition function.
- `ne::Flt`: the number density of free electrons.
- `ν::Flt`: frequency in Hz
- `ρ::Flt`: mass density in g/cm³
- `T::Flt`: temperature in K

# Notes
This is taken from equation 8.13 of Gray (2005). This uses a polynomial to Bell & Berrington (1987)
tabulated values for α_ff(H⁻), which implicitly includes the correction for stimulated emission.
The polynomial fits α_ff(H⁻) in the range 2520 K ≤ T ≤ 10080 K and 2604 Å ≤ λ ≤ 113918 Å. According
to the book, the polynomial fit to the tabulated data typically has 1% precision. We find that at
worst, the discrepancy never exceeds 2.25%.

From equations 8.13, 8.18, and 8.19 of Gray (2005), the free-free opacity from H⁻ is given by:
                              n(H I)          n(H I) + n(H II)
   κ_ν = α_ff(H⁻) * Pₑ * ------------------ * -----------------
                          n(H I) + n(H II)            ρ
This can be rewritten as: κ_ν = α_ff(H⁻) * Pₑ * n(H I) / ρ

Based on Section 5.3 from Kurucz (1970), I'm fairly confident that the "more correct" version of
the opacity equation should actually read  κ_ν = α_ff(H⁻) * Pₑ * n(H I, n = 1) / ρ, and that the
form provided by Gray (2005) implicitly assumes that n(H I, n = 1) ~ n(H I). This seems to be a 
prettygood assumption given the tabulated values of the partition function in Table D.2 of Gray 
(2005). From the boltzmann equation:
     n(H I, n = 1) = n(H I)*gₙ₌₁/u(T)*exp(-Eₙ₌₁/(k*T)) = n(H I) * 2/u(T)*exp(0) = n(H I) * 2/u(T).
According to Table D.2 of Grey, this approximation only introduces a ≤0.5% overestimate in 
n(H I, n = 1) over the temperature range where the polynomial is valid.

We also considered the polynomial fit in Section 5.3 from Kurucz (1970). Unfortunately, it seems
wrong (it gives lots of negative numbers).
"""
function Hminus_ff(nH_I_div_partition::Real, ne::Real, ν::Real, ρ::Real, T::Real)
    λ = c_cgs*1e8/ν # in Angstroms
    # we need to somehow factor out this bounds checking
    if !(2604 <= λ <= 113918.0)
        throw(DomainError(λ, "The wavelength must lie in the interval [2604 Å, 113918 Å]"))
    elseif !(2520.0 <= T <= 10080.0)
        throw(DomainError(T, "The temperature must lie in the interval [2520 K, 10080 K]"))
    end

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

    return αff_H⁻ * Pₑ * nHI_gs / ρ
end
    

"""
    H2plus_bf_and_ff(nH_I_div_partition, n_HII, ν, ρ, T)

Compute the combined H₂⁺ bound-free and free-free opacity κ.

This uses polynomial fits from Gray (2005) that were derived from data tabulated in
[Bates (1952)](https://ui.adsabs.harvard.edu/abs/1952MNRAS.112...40B/abstract).

# Arguments
- `nH_I_div_partition`: the total number density of H I divided by its partition 
   function.
- `nH_II`: the number density of H II (not of H₂⁺).
- `ν`: frequency in Hz
- `ρ`: mass density in g/cm³
- `T`: temperature in K

# Notes
This follows equation 8.15 of Gray (2005), which involves 2 polynomials that were fit to data
provided in [Bates (1952)](https://ui.adsabs.harvard.edu/abs/1952MNRAS.112...40B/abstract).

The combined H₂⁺ bound-free and free-free opacity opacity calculation can be cast as:
    κ = const * σ₁ * exp(-U₁ / (kboltz*T)) * n(H I, n=1) * n(H II) * (1 - exp(-hν/(kboltz*T))) / ρ
where σ₁ and U₁ are just functions of ν. Note that Bates (1952) and Gray (2005) both use the number
density of all energy states of H I instead of n(H I, n=1). However, Kurucz (1970) makes a note in
section 5.2 about how the photodisociation of H₂⁺ produces a ground state H atom and that we should
therefore use n(H I, n=1) instead. This difference will only cause Bates/Gray to be wrong by a
fraction of a percent (at most) for T ≤ 1.2e4 K. In this function we use n(H I, n=1).

In this function, we use polynomial approximations that Gray (2005) fit for σ₁ and U₁ using data
from Table 1 of Bates (1952). Gray (2005) notes that the approximations match Bates' (1952)
tabulated absorption data for 3800 Å ≤ λ ≤ 25000 Å at an accuracy of 0.3%. After reviewing the data
in Bates (1952) I've concluded that Gray (2005) actually rounded down the lower end of the
wavelength from 10⁵/26 Å or ∼3847 Å. (In other words, the approximation is good for
3847 Å ≤ λ ≤ 25000.0 Å).

While this may fit to σ₁ and U₁ to better than 0.3%, comparisons against values from table 2 of
Bates (1952) indicate that the function reproduces the full absorption coeffient at better than
1.5%. Gray (2005) did not provide a temperature range, but Bates (1952) only computed the
absorption coefficient for 2500 K ≤ T ≤ 12000 K. In the text they have a comment that the proton's
De Broglie wavelength is large at 2500 K and their semi-classical treatment may start to break
down. However, as long as use n(H I, n=1), there doesn't seem to be any reason for us to enforce an
upper Temperature limit.

Bates (1952) states that his classical treatment of the interaction is probably most accurate at
higher temperatures and longer wavelengths (note that the longest λ considered in the paper is 8
times larger than the max λ the polynomials are fit against). He suggests that treatment is
probably correct "to well within one part in ten even at the lower temperatures and [lower
wavelengths]."
"""
function H2plus_bf_and_ff(nH_I_div_partition::Real, nH_II::Real, ν::Real, ρ::Real, T::Real)
    λ = c_cgs*1e8/ν # in ångstroms
    if !(3846.15 <= λ <= 25000.0) # the lower limit is set to include 1.e5/26 Å
        throw(DomainError(λ, "The wavelength must lie in the interval [3847 Å, 25000 Å]"))
    end

    if T < 2500  # we might be able to drop upper limit
        throw(DomainError(T, "The temperature must be greater than or equal to 2500 K."))
    end

    β_eV = 1.0/(kboltz_eV * T)

    # coef should be roughly 2.51e-42
    coef = 16*π^4*bohr_radius_cgs^5*electron_charge_cgs^2/(3*hplanck_cgs*c_cgs) # cm⁵
    stimulated_emission_correction = (1 - exp(-hplanck_eV*ν*β_eV))

    logλ = log10(λ)
    log2λ = logλ*logλ
    log3λ = log2λ*logλ

    σ1 = -1040.54 + 1345.71 * logλ - 547.628 * log2λ + 71.9684 * log3λ
    # there was a typo in Gray (2005). In the book they give the following polynomial as the fit
    # for U₁. In reality, it is the fit for negative U₁
    neg_U1 = 54.0532 - 32.713 * logλ + 6.6699 * log2λ - 0.4574 * log3λ
    # note: Gray (2005) used the equivalent expression: σ1*10^(neg_U1 * θ) where θ = 5040/T
    atomic_cross_section = σ1 * exp(neg_U1 * β_eV)

    # see the docstring for an explanation of why we explicitly consider the ground state density
    nH_I_gs = ndens_state_hydrogenic(1, nH_I_div_partition, T, _H_I_ion_energy)

    uncorrected_opacity = coef * atomic_cross_section * nH_I_gs * nH_II/ρ

    # Gray (2005) notes that they remove the stimulated emission factor. We need to put it back:
    uncorrected_opacity * stimulated_emission_correction
end
