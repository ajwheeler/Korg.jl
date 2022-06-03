using Interpolations: LinearInterpolation, Throw
using HDF5: h5read

using ..ContinuumAbsorption: hydrogenic_bf_absorption, hydrogenic_ff_absorption, ionization_energies

const _H_I_ion_energy = ionization_energies[1][1] # not sure if this is a good idea
const _H⁻_ion_energy = 0.754204 # [eV] used by McLaughlin+ 2017 H⁻ ff cross sections

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

const _Hminus_bf_cross_section_interp, _min_H⁻_interp_ν = let
    fn = joinpath(_data_dir, "McLaughlin2017Hminusbf.h5")
    ν = h5read(fn, "nu")
    σ = h5read(fn, "sigma")
    LinearInterpolation(ν, σ, extrapolation_bc=Throw()), minimum(ν)
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
  ``\\alpha_{bf}`` in place of ``\\sigma_{bf}``). We get this from 
   [McLaughlin (2017)](https://ui.adsabs.harvard.edu/abs/2017JPhB...50k4001M/abstract).
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

!!! note
    The McLaughlin tabulated values (downloaded as a text file) contained a stray line out of 
    monotonic order.  A version of the data file with the line deleted is saved at 
    `data/McLaughlin2017Hminusbf.dat` for archival purposes.  (Korg doesn't read this file, it reads
    `data/McLaughlin2017Hminusbf.h5`.)

"""
Hminus_bf = bounds_checked_absorption(
    _Hminus_bf;
    ν_bound = closed_interval(0.0, 2.417989242625068e19),
    temp_bound = Interval(0, Inf)
)

const _Hminus_ff_absorption_interp = let 
    # table from Bell & Berrington (1987) https://doi.org/10.1088/0022-3700/20/4/019 
    # (pinched from MOOG source)
    theta_ff_absorption_interp = [0.5,  0.6, 0.8,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.8,  3.6]
    lambda_ff_absorption_interp = [1823, 2278, 2604, 3038, 3645, 4557, 5063, 5696, 6510, 7595, 9113, 
                                  10126, 11392, 13019, 15189, 18227, 22784, 30378, 45567, 91134,
                                  113918, 151890]
    ff_absorption = [
        .0178 .0222 .0308 .0402 .0498 .0596 .0695 .0795 .0896  .131  .172 
        .0228 .0280 .0388 .0499 .0614 .0732 .0851 .0972  .110  .160  .211
        .0277 .0342 .0476 .0615 .0760 .0908  .105  .121  .136  .199  .262
        .0364 .0447 .0616 .0789 .0966  .114  .132  .150  .169  .243  .318
        .0520 .0633 .0859  .108  .131  .154  .178  .201  .225  .321  .418
        .0791 .0959  .129  .161  .194  .227  .260  .293  .327  .463  .602
        .0965  .117  .157  .195  .234  .272  .311  .351  .390  .549  .711
        .121  .146  .195  .241  .288  .334  .381  .428  .475  .667  .861
        .154  .188  .249  .309  .367  .424  .482  .539  .597  .830  1.07
        .208  .250  .332  .409  .484  .557  .630  .702  .774  1.06  1.36
        .293  .354  .468  .576  .677  .777  .874  .969  1.06  1.45  1.83
        .358  .432  .572  .702  .825  .943  1.06  1.17  1.28  1.73  2.17
        .448  .539  .711  .871  1.02  1.16  1.29  1.43  1.57  2.09  2.60
        .579  .699  .924  1.13  1.33  1.51  1.69  1.86  2.02  2.67  3.31
        .781  .940  1.24  1.52  1.78  2.02  2.26  2.48  2.69  3.52  4.31
        1.11  1.34  1.77  2.17  2.53  2.87  3.20  3.51  3.80  4.92  5.97
        1.73  2.08  2.74  3.37  3.90  4.50  5.01  5.50  5.95  7.59  9.06
        3.04  3.65  4.80  5.86  6.86  7.79  8.67  9.50  10.3  13.2  15.6
        6.79  8.16  10.7  13.1  15.3  17.4  19.4  21.2  23.0  29.5  35.0
        27.0  32.4  42.6  51.9  60.7  68.9  76.8  84.2  91.4  117.  140.
        42.3  50.6  66.4  80.8  94.5  107.  120.  131.  142.  183.  219.
        75.1  90.0  118.  144.  168.  191.  212.  234.  253.  325.  388.
    ]
    LinearInterpolation((lambda_ff_absorption_interp, theta_ff_absorption_interp), ff_absorption; 
                         extrapolation_bc=Throw());
end

function _Hminus_ff(ν::Real, T::Real, nH_I_div_partition::Real, ne::Real)
    λ = c_cgs*1e8/ν # in Angstroms
    θ = 5040.0/T

    # K is the variable used for the quantity in the paper. It's units are cm^4/dyn
    #The factor of 1e-26 is built in to the table.
    K = 1e-26 * _Hminus_ff_absorption_interp(λ, θ) 

    # Pₑ * α_ff(H⁻) gives the absorption coefficient in units of cm² per ground state H I atom
    Pₑ = ne * kboltz_cgs * T

    # Account for the fact that n(H I, n=1) might be slightly smaller than the entire number
    # density of H I. There is only really a difference at the highest temperatures. For the
    # temperature range where this approximation is valid, less than 0.23% of all H I atoms are not
    # in the ground state. This calculation could reduce to nHI_gs = 2.0*nH_I_div_partition
    nHI_gs = ndens_state_hydrogenic(1, nH_I_div_partition, T, _H_I_ion_energy)

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
Since U(T) ≈ 2 up to fairly large tempuratures, this is not unreasonable.
"""
Hminus_ff = bounds_checked_absorption(_Hminus_ff,
                                      ν_bound = λ_to_ν_bound(closed_interval(1823e-8,151890e-8)),
                                      temp_bound = closed_interval(1400, 10080))

include("Stancil1994.jl") #used for H2plus_bf_and_ff
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
