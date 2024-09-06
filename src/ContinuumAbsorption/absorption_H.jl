using Interpolations: linear_interpolation, Throw, Line
using HDF5: h5read, h5open

using ..ContinuumAbsorption: hydrogenic_ff_absorption, ionization_energies

const _H⁻_ion_energy = 0.754204 # [eV] used by McLaughlin+ 2017 H⁻ ff cross sections

# load hydrogen bf cross sections
const _H_I_bf_cross_sections = let
    h5open(joinpath(_data_dir, "bf_cross-sections",
                    "individual_H_cross-sections.h5")) do f
        sigmas = map(eachcol(read(f["E"])), eachcol(read(f["sigma"]))) do Es, σs
            linear_interpolation(Es, σs; extrapolation_bc=Line())
        end
        # use the cross sections for the first 6 energy levels only.
        # the binding energy for n=7 corresponds to ~45,000 Å 
        collect(zip(read(f["n"]), sigmas))
    end
end

"""
    H_I_bf(νs, T, nH, nHe, ne, invU_H; n_max_MHD=6, use_hubeny_generalization=false, 
           taper=false, use_MHD_for_Lyman=false)

The bound-free linear absorption coefficient contributed by all energy states of a neutral Hydrogen
atom. Even though the Mihalas-Hummer-Daeppen (MHD) occupation probability formalism is not used in
Korg when computing the hydrogen partition function, it is used here. That means that series limit
jumps (e.g. the Balmer jump), are "rounded off", as photons with less than the "classical"
ionization energy can ionize if the upper level is dissolved into the continuum.

# Required Arguments

  - `νs`: sorted frequency vector in Hz
  - `T`: temperature in K
  - `nH_I`: the total number density of neutral Hydrogen (in cm⁻³)
  - `nHe_I`: the total number density of neutral Helium (in cm⁻³)
  - `ne`: the number density of electrons (in cm⁻³)
  - `invU_H`: The inverse of the neutral hydrogen partition function (neglecting contributions from the
    MHD formalism)

For n=1 through n=`n_max_MHD` (default: 6), the cross-sections are computed using
[Nahar 2021](https://ui.adsabs.harvard.edu/abs/2021Atoms...9...73N/abstract).
These are modified using the MHD formalism to account for level dissolution.
For larger n, the cross-sections are calculated with a simple analytic formula (see
[`simple_hydrogen_bf_cross_section`](@ref)).

Because MHD level dissolution applied to the the Lyman series limit leads to inflated cross-sections
in the visible, we don't use MHD for bf absorption from n=1.  This can be overridden by setting
`use_MHD_for_Lyman=true`, in which case you will also want to set `taper=true`, which the same
tapering of the cross-section as [HBOP](https://github.com/barklem/hlinop/blob/master/hbop.f) to fix
the problem.

The `use_hubeny_generalization` keyword argument enables the generalization of the MHD from
Hubeny 1994. It is experimental and switched off by default.
"""
function H_I_bf(νs, T, nH_I, nHe_I, ne, invU_H; n_max_MHD=6, use_hubeny_generalization=false,
                taper=false, use_MHD_for_Lyman=false)
    χ = ionization_energies[1][1]
    σ_type = promote_type(eltype(νs), typeof(T), typeof(nH_I), typeof(nHe_I), typeof(ne),
                          typeof(invU_H))
    total_cross_section = zeros(σ_type, length(νs))
    # allocate working vectors outside of loop
    cross_section = Vector{σ_type}(undef, length(νs))
    dissolved_fraction = Vector{σ_type}(undef, length(νs))
    for (n, sigmas) in _H_I_bf_cross_sections[1:n_max_MHD]
        w_lower = hummer_mihalas_w(T, n, nH_I, nHe_I, ne;
                                   use_hubeny_generalization=use_hubeny_generalization)
        #the degeneracy is already factored into the nahar cross-sections
        occupation_prob = w_lower * exp(-χ * (1 - 1 / n^2) / (kboltz_eV * T))

        ν_break = χ / (n^2 * hplanck_eV)
        # the index of the lowest frequency in νs that is greater than the ionization frequency for 
        # an unperturbed H atom with its electron in the nth level   
        break_ind = findfirst(νs .> ν_break)
        if isnothing(break_ind)
            break_ind = length(νs) + 1
        end

        # dissolved_fraction in the fraction of atoms for which the upper state corresponding to ν 
        # is dissolved into the continuum.
        # all the levels above the unperturbed ionization threshold are dissolved by definition
        dissolved_fraction[break_ind:end] .= 1.0
        if !use_MHD_for_Lyman && n == 1
            # don't use MHD for the Lyman series limit since it leads to inflated cross-sections
            # far red of the limit
            dissolved_fraction[1:break_ind-1] .= 0.0
        else
            dissolved_fraction[1:break_ind-1] .= map(νs[1:break_ind-1]) do ν
                # account for bf absorption redward of the jump/break because of level dissolution
                # the effective quantum number associated with the energy of the nth level plus the 
                # photon energy
                n_eff = 1 / sqrt(1 / n^2 - hplanck_eV * ν / χ)
                # this could probably be interpolated without much loss of accuracy
                w_upper = hummer_mihalas_w(T, n_eff, nH_I, nHe_I, ne;
                                           use_hubeny_generalization=use_hubeny_generalization)
                # w_upper/w[n] is the prob that the upper level is dissolved given that the lower isn't
                frac = 1 - w_upper / w_lower
                if taper
                    # taper of the cross-section past a certain wavelength redward of the jump/break 
                    # (non-MHD ionization energy), as is done in HBOP. (Not enabled in Korg calls to 
                    # this function.)
                    redcut = hplanck_eV * c_cgs / (χ * (1 / n^2 - 1 / (n + 1)^2))
                    λ = c_cgs / ν
                    if λ > redcut
                        frac *= exp(-(λ - redcut) * 1e6)
                    end
                end
                frac
            end
        end

        # extrapolate the cross-section to lower energies by assuming proportionality to ν^-3
        σ_break = sigmas(ν_break * hplanck_eV)
        scaling_factor = σ_break / ν_break^-3
        cross_section[1:break_ind-1] .= νs[1:break_ind-1] .^ -3 * scaling_factor
        cross_section[break_ind:end] .= sigmas.(hplanck_eV .* νs[break_ind:end])

        total_cross_section .+= occupation_prob .* cross_section .* dissolved_fraction
    end

    for n in (n_max_MHD+1):40
        w_lower = hummer_mihalas_w(T, n, nH_I, nHe_I, ne;
                                   use_hubeny_generalization=use_hubeny_generalization)
        if w_lower < 1e-5
            break
        end
        occupation_prob = 2n^2 * w_lower * exp(-χ * (1 - 1 / n^2) / (kboltz_eV * T))
        @. total_cross_section += occupation_prob * simple_hydrogen_bf_cross_section.(n, νs)
    end

    #factor of 10^-18 converts cross-sections from megabarns to cm^2
    @. nH_I * invU_H * total_cross_section * (1.0 - exp(-hplanck_eV * νs / (kboltz_eV * T))) * 1e-18
end

"""
    simple_hydrogen_bf_cross_section(n::Integer, ν::Real)

Calculate the H I bf cross section in megabarns using a very simple approximation.  See, for
example, Kurucz 1970 equation 5.5 (though see note below).  This implementation is used to
extrapolate the cross-section past the ionization energy of an unperturbed hydrogen atom, as is
required to take level dissolution into account with the MHD formalism.

Equation 5.5 of Kurucz had a typo in it. In the numerator of the fraction that is multiplied by the
entire polynomial, Z² should be Z⁴. This was discovered during comparisons with data from the
Opacity Project, and can be confirmed by looking at eqn 5.6 of Kurucz (it uses Z⁴ instead of Z²) or
by comparison against equation 10.54 of Rybicki & Lightman.
"""
function simple_hydrogen_bf_cross_section(n::Integer, ν::Real)
    # this implements equation 5.5 from Kurucz (1970)
    # - Z is the atomic number (Z=1 for hydrogen atoms)
    # - n is the energy level (remember, they start from n=1)
    # - ν should have units of Hz
    # - we are using the polynomial coefficients appropriate for n>6
    if (n < 7)
        if (n < 1)
            throw(DomainError(n, "n must be a positive integer"))
        else
            @warn "this function may provide inaccurate estimates when the n < 7"
        end
    end

    inv_n = 1.0 / n
    inv_n2 = inv_n * inv_n
    if hplanck_eV * ν < (RydbergH_eV * inv_n2)
        # photon doesn't carry minimum energy required to ionize the photon
        # -> definitionally the cross-section is zero
        0.0
    else
        inv_ν = 1.0 / ν
        # 64*π⁴e^10 mₑ / (c h⁶ 3√3), where e is the elementary charge in units of statcoulombs
        bf_σ_const = 2.815e29
        bf_σ_const * (inv_n2 * inv_n2 * inv_n) * (inv_ν^3) * 1e18 #convert to megabarns
    end
end

"""
    _ndens_Hminus(nH_I_div_partition, ne, T, ion_energy = _H⁻_ion_energy)

Compute the number density of H⁻ (implements eqn 5.10 of Kurucz 1970). This is an application of
the saha equation where the "ground state" is H⁻ and the "first ionization state" is H I. The
partition function of H⁻ is 1 at all temperatures.
"""
function _ndens_Hminus(nH_I_div_partition, ne, T, ion_energy=_H⁻_ion_energy)
    # the Boltzmann factor is unity for the ground state, and the degeneracy, g, is 2.
    nHI_groundstate = 2 * nH_I_div_partition

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
    β = 1.0 / (kboltz_eV * T)
    0.25 * nHI_groundstate * ne * coef * β^1.5 * exp(ion_energy * β)
end

const _Hminus_bf_cross_section_interp, _min_H⁻_interp_ν = let
    fn = joinpath(_data_dir, "McLaughlin2017Hminusbf.h5")
    ν = h5read(fn, "nu")
    σ = h5read(fn, "sigma")
    linear_interpolation(ν, σ; extrapolation_bc=Throw()), minimum(ν)
end
#returns the cross-section in units of cm² (excludes stimulated emission)
function _Hminus_bf_cross_section(ν)
    _H⁻_ion_ν = _H⁻_ion_energy / hplanck_eV
    if ν <= _H⁻_ion_ν
        0.0
    elseif ν < _min_H⁻_interp_ν
        #McLaughlin+ 2017 notes that for Eᵧ < 0.7678 eV, that they use σ = 460.8*(Eᵧ - E₀)^1.5 Mb, 
        #where E₀ is the ionization energy. This is indeed consistent with the table, but 460.8 is 
        #rounded up. We will use this same scaling, but in terms of frequency.
        _H⁻_low_ν_coef = (_Hminus_bf_cross_section_interp(_min_H⁻_interp_ν) /
                          (_min_H⁻_interp_ν - _H⁻_ion_ν)^1.5)
        _H⁻_low_ν_coef * (ν - _H⁻_ion_ν)^1.5
    else
        _Hminus_bf_cross_section_interp(ν)
    end
end

function _Hminus_bf(ν::Real, T::Real, nH_I_div_partition::Real, ne::Real)
    cross_section = _Hminus_bf_cross_section(ν) # in units of cm²
    stimulated_emission_correction = (1 - exp(-hplanck_cgs * ν / (kboltz_cgs * T)))
    n_Hminus = _ndens_Hminus(nH_I_div_partition, ne, T, _H⁻_ion_energy)
    n_Hminus * cross_section * stimulated_emission_correction
end

"""
    Hminus_bf(ν, T, nH_I_div_partition, ne; kwargs...)

Compute the H⁻ bound-free linear absorption coefficient α,
``\\alpha_\\nu = \\sigma_{bf}(H^-) n(H⁻) (1 - \\exp \\left( \\frac{-h\\nu}{k T}\\right))``

# Arguments

  - `ν::AbstractVector{<:Real}`: sorted frequency vector in Hz
  - `T`: temperature in K
  - `nH_I_div_partition`: the total number density of H I divided by its partition function.
  - `ne`: the electron number density

This uses cross-sections from
[McLaughlin 2017](https://ui.adsabs.harvard.edu/abs/2017JPhB...50k4001M/abstract).
For a description of the kwargs, see [Continuum Absorption Kwargs](@ref).

# Notes

This function assumes that n(H⁻) ≪ n(H I) + n(H II). The number density of n(H⁻) isn't precomputed
as part of Korg's molecular equlibrium, it's computed here instead.

!!! note

    The McLaughlin tabulated values (downloaded as a text file) contained a stray line out of
    monotonic order.  A version of the data file with the line deleted is saved at
    `data/McLaughlin2017Hminusbf.dat` for archival purposes.  (Korg doesn't read this file, it reads
    `data/McLaughlin2017Hminusbf.h5`.)
"""
Hminus_bf = bounds_checked_absorption(_Hminus_bf;
                                      ν_bound=closed_interval(0.0, 2.417989242625068e19),
                                      temp_bound=Interval(0, Inf))

const _Hminus_ff_absorption_interp = let
    # table from Bell & Berrington (1987) https://doi.org/10.1088/0022-3700/20/4/019 
    theta_ff_absorption_interp = [0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.8, 3.6]
    lambda_ff_absorption_interp = [1823, 2278, 2604, 3038, 3645, 4557, 5063, 5696, 6510, 7595, 9113,
        10126, 11392, 13019, 15189, 18227, 22784, 30378, 45567, 91134,
        113918, 151890]
    ff_absorption = [0.0178 0.0222 0.0308 0.0402 0.0498 0.0596 0.0695 0.0795 0.0896 0.131 0.172
                     0.0228 0.0280 0.0388 0.0499 0.0614 0.0732 0.0851 0.0972 0.110 0.160 0.211
                     0.0277 0.0342 0.0476 0.0615 0.0760 0.0908 0.105 0.121 0.136 0.199 0.262
                     0.0364 0.0447 0.0616 0.0789 0.0966 0.114 0.132 0.150 0.169 0.243 0.318
                     0.0520 0.0633 0.0859 0.108 0.131 0.154 0.178 0.201 0.225 0.321 0.418
                     0.0791 0.0959 0.129 0.161 0.194 0.227 0.260 0.293 0.327 0.463 0.602
                     0.0965 0.117 0.157 0.195 0.234 0.272 0.311 0.351 0.390 0.549 0.711
                     0.121 0.146 0.195 0.241 0.288 0.334 0.381 0.428 0.475 0.667 0.861
                     0.154 0.188 0.249 0.309 0.367 0.424 0.482 0.539 0.597 0.830 1.07
                     0.208 0.250 0.332 0.409 0.484 0.557 0.630 0.702 0.774 1.06 1.36
                     0.293 0.354 0.468 0.576 0.677 0.777 0.874 0.969 1.06 1.45 1.83
                     0.358 0.432 0.572 0.702 0.825 0.943 1.06 1.17 1.28 1.73 2.17
                     0.448 0.539 0.711 0.871 1.02 1.16 1.29 1.43 1.57 2.09 2.60
                     0.579 0.699 0.924 1.13 1.33 1.51 1.69 1.86 2.02 2.67 3.31
                     0.781 0.940 1.24 1.52 1.78 2.02 2.26 2.48 2.69 3.52 4.31
                     1.11 1.34 1.77 2.17 2.53 2.87 3.20 3.51 3.80 4.92 5.97
                     1.73 2.08 2.74 3.37 3.90 4.50 5.01 5.50 5.95 7.59 9.06
                     3.04 3.65 4.80 5.86 6.86 7.79 8.67 9.50 10.3 13.2 15.6
                     6.79 8.16 10.7 13.1 15.3 17.4 19.4 21.2 23.0 29.5 35.0
                     27.0 32.4 42.6 51.9 60.7 68.9 76.8 84.2 91.4 117.0 140.0
                     42.3 50.6 66.4 80.8 94.5 107.0 120.0 131.0 142.0 183.0 219.0
                     75.1 90.0 118.0 144.0 168.0 191.0 212.0 234.0 253.0 325.0 388.0]
    linear_interpolation((lambda_ff_absorption_interp, theta_ff_absorption_interp), ff_absorption;
                         extrapolation_bc=Throw())
end

function _Hminus_ff(ν::Real, T::Real, nH_I_div_partition::Real, ne::Real)
    λ = c_cgs * 1e8 / ν # in Angstroms
    θ = 5040.0 / T

    # K is the variable used for the quantity in the paper. It's units are cm^4/dyn
    #The factor of 1e-26 is built in to the table.
    K = 1e-26 * _Hminus_ff_absorption_interp(λ, θ)

    # Pₑ * α_ff(H⁻) gives the absorption coefficient in units of cm² per ground state H I atom
    Pₑ = ne * kboltz_cgs * T

    # Account for the fact that n(H I, n=1) might be slightly smaller than the entire number
    # density of H I. There is only really a difference at the highest temperatures. For the
    # temperature range where this approximation is valid, less than 0.23% of all H I atoms are not
    # in the ground state. The Boltzmann factor is unity and the degeneracy is 2 for the ground 
    # state.
    nHI_gs = 2 * nH_I_div_partition

    return K * Pₑ * nHI_gs
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

This is based on Table 1, in
[Bell & Berrington (1987)](https://doi.org/10.1088/0022-3700/20/4/019), which tabulates values
for "the H⁻ absorption coefficient", K (including the correction for stimulated emission). This
quantity is in units of cm^4/dyn, and must be multiplied by the electron partial pressure and the
ground-state neutral hydrogen number density to obtain a linear absorption coefficent, α.

The stipulation that the hydrogen should be ground-state only is based on the beginning of Section 2
in Bell and Berrington (1987), or alternately, Section 5.3 from Kurucz (1970).  When Gray (2005)
refers to this, it implicitly assumes that `n(H I, n = 1) ≈ n(H I)`.  Note that

```
 n(H I, n = 1) = n(H I)*gₙ₌₁/U(T)*exp(-Eₙ₌₁/(k*T)) = n(H I) * 2/U(T)*exp(0) = n(H I) * 2/U(T).
```

Since U(T) ≈ 2 up to fairly large temperatures, this is not unreasonable.
"""
Hminus_ff = bounds_checked_absorption(_Hminus_ff;
                                      ν_bound=λ_to_ν_bound(closed_interval(1823e-8, 151890e-8)),
                                      temp_bound=closed_interval(1400, 10080))

include("Stancil1994.jl") #used for H2plus_bf_and_ff
function _H2plus_bf_and_ff(ν::Real, T::Real, nH_I::Real, nH_II::Real)
    λ = c_cgs * 1e8 / ν # in ångstroms
    β_eV = 1.0 / (kboltz_eV * T)
    stimulated_emission_correction = (1 - exp(-hplanck_eV * ν * β_eV))

    K = Stancil1994.K_H2plus(T) # n(H I) n(H II) / n(H₂⁺)
    σbf = Stancil1994.σ_H2plus_bf(λ, T)
    σff = Stancil1994.σ_H2plus_ff(λ, T)

    (σbf / K + σff) * nH_I * nH_II * stimulated_emission_correction
end

"""
    H2plus_bf_and_ff(ν, T, nH_I_div_partition, n_HII; kwargs...)

Compute the combined H₂⁺ bound-free and free-free linear absorption coefficient α using the tables
from [Stancil 1994](https://ui.adsabs.harvard.edu/abs/1994ApJ...430..360S/abstract). Note that
this refers to interactions between free-protons and neutral hydrogen, not electrons and doubly
ionized H2, i.e. the b-f interaction is photodissociation, not photoionization.

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
H2plus_bf_and_ff = bounds_checked_absorption(_H2plus_bf_and_ff;
                                             ν_bound=λ_to_ν_bound(closed_interval(7e-6, 2e-3)),
                                             temp_bound=closed_interval(3150, 25200))
