using Interpolations: LinearInterpolation, Throw
using StaticArrays: SA

# define the hydrogenic continuum absorption

"""
    ndens_state_hydrogenic(n, nsdens_div_partition, ion_energy, T)

Calculates the LTE number density (in cm^-3) of a hydrogenic species at a given energy level.

# Arguments
- `n::Integer`: The quantum number of the state
- `nsdens_div_partition`: The total (including all states) number density of the current 
ionization species (in cm^-3) divided by the species's partition function.
- `ion_energy`: The minimum energy (in eV) required to ionize the species (from the ground 
state). This can be estimated as Z²*ion_energy of hydrogen or Z²*Rydberg_H.
- `T`: Temperature

# Note
This assumes that all hydrogenic species have a statistical weight of gₙ = 2*n².

This was taken from equation (5.4) of Kurucz (1970) (although this comes directly from the 
boltzmann equation).
"""
function ndens_state_hydrogenic(n::Integer, nsdens_div_partition::Real, T::Real, ion_energy::Real)
    n2 = n*n
    g_n = 2.0*n2
    energy_level = ion_energy - ion_energy/n2
    nsdens_div_partition * g_n * exp(-energy_level/(kboltz_eV *T))
end

_eVtoHz(energy) = energy/hplanck_eV
_HztoeV(freq) = freq*hplanck_eV

# In the following equation, e is the elementary charge in units of statcoulombs
const _bf_σ_const = 2.815e29 # = 64*π⁴e^10 mₑ / (c h⁶ 3√3)

const _bf_σ_coef = SA[(0.9916,  2.719e13, -2.268e30),
                      (1.105,  -2.375e14,  4.077e28),
                      (1.101,  -9.863e13,  1.035e28),
                      (1.101,  -5.765e13,  4.593e27),
                      (1.102,  -3.909e13,  2.371e27),
                      (1.0986, -2.704e13,  1.229e27),
                      (1.0,     0.0,       0.0     )]
 
# this is helper function is the most likely part of the calculation to change.
# This uses double precision to be safe about the polynomial coefficients.
function _hydrogenic_bf_cross_section(Z::Integer, n::Integer, ν::Real, ion_freq::Real)
    # this implements equation 5.5 from Kurucz (1970)
    # - Z is the atomic number
    # - n is the energy level (remember, they start from n=1)
    # - ν should have units of Hz
    # - ion_freq is the frequency of the photon carrying the minimum energy needed to ionize the
    #   ground state configuration of the current ion (in Hz). This can be estimated as
    #   _eVtoHz(Z²*RydbergH_eV).

   if (n < 1)
        throw(DomainError(n,"n must be a positive integer"))
    end

    Z2 = Z * Z
    inv_n = 1.0/n
    inv_n2 = inv_n*inv_n
    inv_ν = 1.0/ν
    if ν < (ion_freq*inv_n2)
        0.0
    else
        # the last entry of the table should be used when n > length(_bf_σ_coef)
        A, B, C = _bf_σ_coef[min(n,length(_bf_σ_coef))]

        poly = A + (B + C*Z2*inv_ν)*Z2*inv_ν
        # equation 5.5 of Kurucz had a typo in it. In the numerator of the fraction that is
        # multiplied by the entire polynomial, Z² should be Z⁴. This was discovered during
        # comparisons with data from the Opacity Project, and can be confirmed by:
        # - looking at eqn 5.6 of Kurucz (it uses Z⁴ instead of Z²)
        # - comparison against equation 10.54 of Rybicki & Lightman
        _bf_σ_const*Z2*Z2*(inv_n2*inv_n2*inv_n)*(inv_ν^3)*poly
    end
end

"""
    _hydrogenic_bf_high_n_absorption(ν, T, Z, nmin, nsdens_div_partition, ion_energy)

Approximates the sum of absorption coefficient contributions by all bound-free transitions in which
the electron originates in an energy level of nmin or larger using an integral.

In full detail, the calculated quantity is described as
    ∑_{n'=nmin}^∞ ndens(n=n') * α_ν(n=n') * (1 - exp(-hν/k/T))/ρ,
where α(n=n') is the absorption coefficient for the bound-free atomic absorption coefficient
(uncorrected for stimulated emission). This function approximates this sum with an integral.

# Arguments
- `ν`: frequency in Hz
- `T`: temperature in K
- `Z::Integer`: Z is the atomic number of the ion (1 for HI)
- `nmin::Integer`: The lowest energy level (principle quantum number) included in the calculation
- `nsdens_div_partition` is the number density of the current species divided by the
   partition function.
- `ion_energy`: the ionization energy from the ground state (in eV).

# Notes
This is based on equation (5.6) from Kurucz (1970) (which give the formula for opacity). I think ρ
was simply omitted from that equation.
"""
function _hydrogenic_bf_high_n_absorption(ν::Real, T::Real, Z::Integer, nmin::Integer,
                                          nsdens_div_partition::Real, ion_energy::Real)

    # this function corresponds to the evaluation of a integral. We subdivide the solution into two
    # parts: (i) consts and (ii) integral. The solution is the product of both parts

    # first, compute consts (terms pulled out of the integral)
    β = 1.0/(kboltz_eV*T)
    expν_term = exp(-β*hplanck_eV*ν)
    integral_consts = nsdens_div_partition * _bf_σ_const * Z^4 * (1.0 - expν_term) * 2.0 / ν^3

    # second, compute the integral
    # need to know the ionization frequency (frequency of a photon that could ionize species
    # from ground state)
    ion_freq = _eVtoHz(ion_energy)

    ionE_times_β = ion_energy*β
    coef = 0.5/ionE_times_β

    # take 2 different approaches based on the photon's energy
    if ν >= (ion_freq/nmin^2)
        # At this frequency, photons are energetic enough to ionize electrons at n >= nmin
        # integral of exp(-ion_energy*(1-(1/n²)) /(kboltz*T))/n³ from nmin to ∞
        integral = coef * (exp(-(ionE_times_β - ionE_times_β/nmin^2)) - exp(-ionE_times_β))
    else
        # At this frequency, photon aren't energetic enough to ionize electrons at:
        # nmin ≤ n < ceil(sqrt(ion_freq/ν)). The absorption coefficients for these electrons are
        # all zero.
        # Thus, the integral's lower limit should actually be ceil(sqrt(ion_freq/ν)). Kurucz
        # approximates this by replacing ion_energy*β/nmin² in the prior solution with hplanck*ν*β
        integral = coef * (exp(-(ionE_times_β - hplanck_eV*ν*β)) - exp(-ionE_times_β))
    end

    (integral_consts*integral)
end


"""
    hydrogenic_bf_absorption(ν, T, Z, nmax_explicit_sum, nsdens_div_partition,
                             ion_energy, [integrate_high_n])

Compute the bound-free linear absorption coefficient contributed by all energy states of a
Hydrogenic species

The calculation is broken into 2 parts: (i) the contributions of the lowest energy states are 
explicitly summed and (ii) the contributions of the higher energy states are estimated with an 
integral.

# Arguments
- `ν`: frequency in Hz
- `T`: temperature in K
- `Z::Integer`: Z is the atomic number of the species (e.g. 1 for H I or 2 for He II)
- `nsdens_div_partition` is the total number density of the species divided by the species's
   partition function.
- `ion_energy`: the ionization energy from the ground state (in eV). This can be 
   estimated as Z²*Rydberg_H (Rydberg_H is the ionization energy of Hydrogen)
- `nmax_explicit_sum::Integer`: The highest energy level whose absorption contribution is included
   in the explicit sum. The contributions from higher levels are included in the integral.
- `integrate_high_n::bool`: When this is `false`, bf absorption from higher energy states are not
   estimated at all. Default is `true`.

# Notes
This follows the approach described in section 5.1 of Kurucz (1970).
"""
function hydrogenic_bf_absorption(ν::Real, T::Real, Z::Integer, nsdens_div_partition::Real,
                                  ion_energy::Real, nmax_explicit_sum::Integer, 
                                  integrate_high_n::Bool=true)
    ionization_freq = _eVtoHz(ion_energy)

    # first, directly sum the individual absorption contributions from H I atoms at each of the
    # lowest energy levels (i.e. all energy levels where n <= nmax_explicit_sum)
    partial_sum = 0.0
    for n = 1 : nmax_explicit_sum
        ndens_state = ndens_state_hydrogenic(n, nsdens_div_partition, T, ion_energy)
        hydrogenic_bf_cross_section = _hydrogenic_bf_cross_section(Z, n, ν, ionization_freq)
        cur_val = ndens_state * hydrogenic_bf_cross_section
        partial_sum += ndens_state * hydrogenic_bf_cross_section
    end
    α_low_n = partial_sum * (1.0 - exp(-hplanck_eV * ν / (kboltz_eV * T)))

    # second, estimate the absorption contributions from H I atoms at higher energy levels using an
    # integral approximation (assuming integrate_high_n is true)
    α_high_n = _hydrogenic_bf_high_n_absorption(ν, T, Z, nmax_explicit_sum+1, nsdens_div_partition,
                                                ion_energy)

    α_low_n + (α_high_n * integrate_high_n)
end



# this table is taken from section 5.1 of Kurucz (1970)
const _ff_interpolator = begin
    table_val = [5.53 5.49 5.46 5.43 5.40 5.25 5.00 4.69 4.48 4.16 3.85;
                 4.91 4.87 4.84 4.80 4.77 4.63 4.40 4.13 3.87 3.52 3.27;
                 4.29 4.25 4.22 4.18 4.15 4.02 3.80 3.57 3.27 2.98 2.70;
                 3.64 3.61 3.59 3.56 3.54 3.41 3.22 2.97 2.70 2.45 2.20;
                 3.00 2.98 2.97 2.95 2.94 2.81 2.65 2.44 2.21 2.01 1.81;
                 2.41 2.41 2.41 2.41 2.41 2.32 2.19 2.02 1.84 1.67 1.50;
                 1.87 1.89 1.91 1.93 1.95 1.90 1.80 1.68 1.52 1.41 1.30;
                 1.33 1.39 1.44 1.49 1.55 1.56 1.51 1.42 1.33 1.25 1.17;
                 0.90 0.95 1.00 1.08 1.17 1.30 1.32 1.30 1.20 1.15 1.11;
                 0.45 0.48 0.52 0.60 0.75 0.91 1.15 1.18 1.15 1.11 1.08;
                 0.33 0.36 0.39 0.46 0.59 0.76 0.97 1.09 1.13 1.10 1.08;
                 0.19 0.21 0.24 0.28 0.38 0.53 0.76 0.96 1.08 1.09 1.09]
    # The x-axis is log₁₀(RydbergH*Z²/(k*T)) = log₁₀(γ²)
    xaxis = [-3.0, -2.5, -2.0, -1.5, -1.0, -0.5,  0.0,  0.5,  1.0,  1.5,  2.0]
    # The y-axis is log₁₀(h*ν/(k*T)) = log₁₀(u)
    yaxis = [-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5,  0.0,  0.5,  1.0,  1.5]

    # Since the grid is regularly spaced, we could probably make this faster
    # pass y,x to the interpolator
    LinearInterpolation((yaxis, xaxis), table_val, extrapolation_bc=Throw())
end

"""
    gaunt_ff_kurucz(log_u, log_γ2)

computes the thermally averaged free-free gaunt factor by interpolating the table provided in
section 5.1 of Kurucz (1970). The table was derived from a figure in Karsas and Latter (1961).

# Arguments
- `log_u`: Equal to log₁₀(u) = log₁₀(h*ν/(k*T))
- `log_γ2`: Equal to log₁₀(γ²) = log₁₀(RydbergH*Z²/(k*T))

# Note
There is some ambiguity over whether we should replace RydbergH*Z² in the definition of γ² with the
ionization energy, but it's probably a negligible difference (given the log scale).

In the future, we may want to the interpolate the table reported in van Hoof, Ferland, Williams,
Volk, Chatzikos, Lykins, & Porter (2013) because it's more accurate and applies over a larger range
of values.
"""
gaunt_ff_kurucz(log_u, log_γ2) = _ff_interpolator(log_u, log_γ2)


"""
    hydrogenic_ff_absorption(ν, T, Z, ni, ne)

computes the free-free linear absorption coefficient for a hydrogenic species

The naming convention for free-free absorption is counter-intuitive. A free-free interaction is
named as though the species interacting with the free electron had one more bound electron (in 
other words it's named as though the free-electron and ion were bound together). In practice, this
means that `ni` should refer to:
- the number density of H II if computing the H I free-free absorption
- the number density of He III if computing the He II free-free absorption
- the number density of Li IV if computing the Li III free-free absorption

# Arguments
- `Z::Integer`: the charge of the ion. For example, this is 1 for ionized H.
- `ni`: the number density of the ion species in cm⁻³.
- `ne`: the number density of free electrons.
- `ν`: frequency in Hz
- `T`: temperature in K

# Note
This approach was adopted from equation 5.8 from section 5.1 of Kurucz (1970). Comparison against
equation 5.18b of Rybicki & Lightman (2004), reveals that the equation in Kurucz (1970) omits the
dependence on ρ. According to Rybicki & Lightman (2004) the free-free absorption coefficient
(corrected for stimulated emission) is:
```
    α = coef * Z² *ne * ni * (1 - exp(-hplanck*ν/(kboltz*T))) * g_ff / (sqrt(T) * ν³)
```
Note that the g_ff is the free-free gaunt factor and coef is ∼3.7e8 (a more exact coefficient can
be computed from eqn 5.18a).

With this in mind, equation 5.8 of Kurucz (1970) should actually read
```
    κ = ne * n(H II) * F_ν(T) * (1 - exp(-hplanck*ν/(kboltz*T))) / ρ
```
where F_ν(T) = coef * Z² * g_ff / (sqrt(T) * ν³).
"""
function hydrogenic_ff_absorption(ν::Real, T::Real, Z::Integer, ni::Real, ne::Real)
    inv_T = 1.0/T
    Z2 = Z*Z

    hν_div_kT = (hplanck_eV/kboltz_eV) * ν * inv_T
    log_u = log10(hν_div_kT)
    log_γ2 = log10((RydbergH_eV/kboltz_eV) * Z2 * inv_T)
    gaunt_ff = gaunt_ff_kurucz(log_u, log_γ2)

    F_ν = 3.6919e8*gaunt_ff*Z2*sqrt(inv_T)/(ν*ν*ν)

    ni*ne*F_ν*(1-exp(-hν_div_kT))
end
