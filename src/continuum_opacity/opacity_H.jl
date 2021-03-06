# we are missing Rayleigh Scattering

using ..ContinuumOpacity: hydrogenic_bf_opacity, hydrogenic_ff_opacity, ionization_energies

const _H_I_ion_energy = ionization_energies["H"][1] # not sure if this is a good idea
const _H⁻_binding_energy = 0.7552 # eV


H_I_bf_opacity(nH_I_div_partition, ν, ρ, T, ion_energy = _H_I_ion_energy) =
               hydrogenic_bf_opacity(1, 8, nH_I_div_partition, ν, ρ, T, ion_energy)
H_I_ff_opacity(nH_I, ne, ν, T) = hydrogenic_ff_opacity(1, nH_I, ne, ν, T)

# compute the number density of H⁻ (implements eqn 5.10 of Kurucz 1970). This formula comes from
# inverting the saha equation, where n(H⁻) is n₀ and n(H I) is n₁. Note that U₀ = 1 at all
# temperatures.
function ndens_Hminus(nH_I_div_partition, ne, T, binding_energy = _H⁻_binding_energy)
    nHI_groundstate = ndens_state_hydrogenic(1, nH_I_div_partition, T, _H_I_ion_energy)

    # coef = (h^2/(2*π*m))^1.5
    coef = 3.31283018e-22 # cm³*eV^1.5
    β = 1.0/(kboltz_eV*T)

    0.25 * nHI_groundstate * ne * coef * β^1.5 * exp(binding_energy * β)
end

function _Hminus_bf_cross_section(ν)
    # taken from Kurucz (1970), who took it from Gingerich (1964)
    inv_ν = 1.0/ν

    if ν >= 2.111e14
        6.801e-20 + (5.358e-3 + (1.481e13 + (-5.519e27 + 4.808e41*inv_ν)*inv_ν)*inv_ν)*inv_ν
    elseif 2.111e14 > ν >= 1.8259e14
        3.695e-16 + (-1.251e-1 + 1.052e13*inv_ν)*inv_ν
    else
        0.
    end
end

"""
    Hminus_bf(nHminus, ν, ρ, T, [binding_energy])

Compute the H⁻ bound-free opacity κ

# Arguments
- `nH_I_div_partition::Flt`: the total number density of H I divided by its partition function.
- `ν::Flt`: frequency in Hz
- `ρ::Flt`: mass density in g/cm³
- `T::Flt`: temperature in K
- `binding_energy::Flt`: Specifies the binding energy of the single state of H⁻ in eV. This is 
   roughly 0.7552 eV.

# Notes
This is taken from Section 5.3 of Kurucz (1970).
The number density of Hminus comes from eqn (5.10) of Kurucz (1970).

An alternative approach is presented in equation 8.11 of Grey (2005) who fit a polynomial to more 
recent work.
"""
function Hminus_bf(nH_I_div_partition::Flt, ne::Flt, ν::Flt, ρ::Flt, T::Flt,
                   binding_energy::Flt = 0.7552) where {Flt<:AbstractFloat}
    nHminus = ndens_Hminus(nH_I_div_partition, ne, T, binding_energy, bHminus)
    one_minus_exp = 1 - exp(-hplanck_cgs*ν/(kboltz_cgs*T))
    α_ν = _Hminus_boundfree_cross_section(ν) # we use alpha here to match notation from the paper
    nHminus * α_ν * one_minus_exp/ρ
end

"""
    Hminus_ff(nHI_gs, ne, ν, ρ, T)

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
This can be rewritten as: κ_ν = α_ff(H⁻) * Pₑ * n(H I) / ρ. We can rearrange this to get:
α_ff(H⁻) = κ_ν * ρ / (n(H I) * Pₑ).

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
function Hminus_ff(nH_I_div_partition::Flt, ne::Flt, ν::Flt, ρ::Flt,
                   T::Flt) where {Flt<:AbstractFloat}
    λ = clight_cgs*1e8/ν # in Angstroms
    # we need to somehow factor out this bounds checking
    if (λ < 2604.0) || (λ > 113918.0)
        throw(DomainError(ν, "The wavelength must lie in the interval [2604 Å, 113918 Å]"))
    elseif (T < 2520.0) || (T > 10080.0)
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

    # in this case, I'm going to account for the fact that n(H I, n=1) might be slightly
    # smaller than the entire number density of H I. This is only a difference at the highest
    # temperatures. For the temperature range where this approximation is valid, less that
    # 0.23% of all H I atoms are not in the ground state.

    # this calculation could reduce to nHI_gs = 2.0*nH_I_div_partition
    nHI_gs = ndens_state_hydrogenic(1, nH_I_div_partition, T, _H_I_ion_energy)

    return αff_H⁻ * Pₑ * nHI_gs / ρ
end
    


"""
    H2plus_bf_and_ff(nH_I_div_partition, nHII, ν, ρ, T)

Compute the combined H₂⁺ bound-free and free-free opacity κ

# Arguments
- `nH_I_div_partition::Float64`: the total number density of H I divided by its partition 
   function.
- `nHII::Float64`: the number density of HII (not of H₂⁺).
- `ν::Float64`: frequency in Hz
- `ρ::Float64`: mass density in g/cm³
- `T::Float64`: temperature in K

# Notes
This is taken from Section 5.2 of Kurucz (1970). They have a whole discussion about the fact that 
since n(H₂⁺) << n(HI) + n(HII) + n(H₂), we can compute the H₂⁺ opacity from n(HI) and n(HII)

This needs to use double precision floats because of the polynomial coefficients
"""
function H2plus_bf_and_ff(nH_I_div_partition::Float64, nHII::Float64, ν::Float64, ρ::Float64,
                          T::Float64)

    # compute the number density of H I in the ground state
    nHI_gs = ndens_state_hydrogenic(1, nH_I_div_partition, T, _H_I_ion_energy)

    one_minus_exp = 1 - exp(-hplanck_cgs*ν/(kboltz_cgs*T))
    # E_s is in units of eV
    E_s = -7.342e-3 + ν* (-2.409e15 +
                          ν * (1.028e-30 +
                               ν * (1.028e-30 +
                                    ν * (-4.230e-46 +
                                         ν * (1.224e-61 - 1.351e-77 * ν)))))
    lnν = log(ν)
    tmp = (-E_s/(kboltz_eV*T))
    F_ν = exp(tmp - 3.0233e3 + lnν*(3.7797e2 + lnν*(-1.82496e1 +
                                                    lnν*(3.9207e-1 - 3.1672e-3 * lnν))))
    nHI_gs * nHII * F_ν * one_minus_exp / ρ
end
