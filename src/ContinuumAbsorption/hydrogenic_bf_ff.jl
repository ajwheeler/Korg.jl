using Interpolations: LinearInterpolation, Throw
using StaticArrays: SA

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
        # calculate the LTE number density of H I atoms in the n=n state
        g_n = 2.0*n^2
        energy_level = ion_energy - ion_energy/n^2
        ndens_state = nsdens_div_partition * g_n * exp(-energy_level/(kboltz_eV *T))

        hydrogenic_bf_cross_section = _hydrogenic_bf_cross_section(Z, n, ν, ionization_freq)
        partial_sum += ndens_state * hydrogenic_bf_cross_section
    end
    α_low_n = partial_sum * (1.0 - exp(-hplanck_eV * ν / (kboltz_eV * T)))

    # second, estimate the absorption contributions from H I atoms at higher energy levels using an
    # integral approximation (assuming integrate_high_n is true)
    α_high_n = _hydrogenic_bf_high_n_absorption(ν, T, Z, nmax_explicit_sum+1, nsdens_div_partition,
                                                ion_energy)

    α_low_n + (α_high_n * integrate_high_n)
end


"""
    _load_gauntff_table([fname])

Returns a table of thermally-averaged free-free Gaunt factors, and the values of log₁₀(γ²) and
log₁₀(u) associated with each point.

This loads the non-relativistic free-free data published by
[van Hoof et al. (2014)](https://ui.adsabs.harvard.edu/abs/2014MNRAS.444..420V/abstract).

Note: This function code could trivially be adapted to load the relativistic free-free gaunt
factors published by
[van Hoof et al (2015)](https://ui.adsabs.harvard.edu/abs/2015MNRAS.449.2112V/abstract).
"""
function _load_gauntff_table(fname = joinpath(_data_dir, "vanHoof2014-nr-gauntff.dat"))
    parse_header_line(type, s) = (last = something(findfirst('#',s), length(s)+1) - 1;
                                  parse.(type,split(s[1:last])) )

    f = open(fname, "r")
    skipchars((c)->false, f; linecomment='#') # skip 1st section of comments

    log10_γ2, log10_u = begin # parse header values
        magic_number = parse_header_line(Int64, readline(f))[1]
        @assert magic_number == 20140210
        num_γ2, num_u = parse_header_line(Int64, readline(f))
        log10_γ2_start = parse_header_line(Float64, readline(f))[1]
        log10_u_start = parse_header_line(Float64, readline(f))[1]
        # step_size has units of dex and applies to both axes
        step_size = parse_header_line(Float64, readline(f))[1]

        (range(log10_γ2_start, length = num_γ2, step = step_size),
         range(log10_u_start, length = num_u, step = step_size))
    end
    skipchars((c)->false, f; linecomment='#') # skip 2nd section of comments

    out = Matrix{Float64}(undef, length(log10_u), length(log10_γ2))
    for i in 1:length(log10_u)
        view(out,i,:) .= parse.(Float64, split(readline(f)))
    end

    # at this point there is a number section of comments followed by a grid of uncertainties for
    # each value. We'll ignore those. That secondary section may or may not be present in versions
    # of the table that included relativistic calculations
    close(f)

    out, log10_γ2, log10_u
end


# maybe also initialize λ_bounds and T_bounds
const _gauntff_interpolator, _gauntff_T_bounds, _gauntff_λ_bounds = begin

    # load in the tabulated data
    table_val, log10_γ2, log10_u = _load_gauntff_table()

    # the table is bigger than it needs to be. We can easily cut its size by a factor of ~10
    T_extrema = [100.0, 1e6] # units of K
    λ_extrema = [1.0e-6, 1.0e-2] # units of cm (equivalent to 100 Å and 100 μm)
    Z_extrema = [1, 2]

    function _find_bound_inds(range, min, max) # this is very similar to some duplicated code
        lb,ub = searchsortedlast(range, min), searchsortedfirst(range, max)
        @assert (range[lb] ≤ min) && (range[ub] ≥ max)
        lb, ub
    end
    calc_log10_γ2(Z, Tₑ) = log10(Rydberg_eV * Z^2 / (kboltz_eV * Tₑ))
    calc_log10_u(λ, Tₑ) = log10(hplanck_cgs * c_cgs / (λ * kboltz_cgs * Tₑ))

    γ2_lb, γ2_ub = _find_bound_inds(log10_γ2, extrema(calc_log10_γ2.(Z_extrema, T_extrema'))...)
    u_lb, u_ub   = _find_bound_inds(log10_u,  extrema( calc_log10_u.(λ_extrema, T_extrema'))...)

    # make a copy of the selected data so that the unneeded data can be garbage collected
    revised_table = copy(table_val[u_lb:u_ub, γ2_lb:γ2_ub])

    (LinearInterpolation((log10_u[u_lb:u_ub], log10_γ2[γ2_lb:γ2_ub]), revised_table,
                         extrapolation_bc = Throw()),
     closed_interval(T_extrema...),
     closed_interval(λ_extrema...))
end

"""
    gaunt_ff_vanHoof(log_u, log_γ2)

computes the thermally averaged, non-relativistic free-free gaunt factor by interpolating the table
provided by [van Hoof et al. (2014)](https://ui.adsabs.harvard.edu/abs/2014MNRAS.444..420V).

# Arguments
- `log_u`: Equal to log₁₀(u) = log₁₀(h*ν/(k*Tₑ))
- `log_γ2`: Equal to log₁₀(γ²) = log₁₀(Rydberg*Z²/(k*Tₑ))

`Rydberg` is the "infinite mass unit of energy" and `Tₑ` is the temperature of free electrons (for
our purposes, we assume that free electrons are in thermal equilibrium with ions and neutral
species).

# Notes

van Hoof et al. (2014) computed the associated data table with a non-relativistic approach, which
is invalid at very high temperatures. They conclude (from comparisons with a different paper) that
their "results should be accurate up to electron temperatures of roughly 100 MK". This is more than
adequate for stellar atmospheres.
In [van Hoof et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015MNRAS.449.2112V), they find that
relativistic effects introduce a ∼0.75% at 100MK, for Z = 1 (when Z > 1, the change is smaller).

This function currently uses linear interpolation. However, van Hoof et al. (2014) provides an
implementation of a third-order Lagrange scheme, which "reaches a relative precision better than
1.5e-4 everywhere." The C and Fortran implementations of this scheme can be found
[here](http://data.nublado.org/gauntff/), and are copyrighted by a BSD-style license.

Earlier variants of this function used less-accurate data from section 5.1 of Kurucz (1970) that
extended over a smaller interval of data. That table was originally derived from a figure in Karsas
and Latter (1961) and it's now used for testing purposes.
"""
gaunt_ff_vanHoof(log_u, log_γ2) = _gauntff_interpolator(log_u, log_γ2)

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
    α = coef * Z² * ne * ni * (1 - exp(-hplanck*ν/(kboltz*T))) * g_ff / (sqrt(T) * ν³)
```
Note that the g_ff is the free-free gaunt factor and coef is ∼3.7e8 (a more exact coefficient can
be computed from eqn 5.18a).

With this in mind, equation 5.8 of Kurucz (1970) should actually read
```
    κ = ne * n(H II) * F_ν(T) * (1 - exp(-hplanck*ν/(kboltz*T))) / ρ
```
where F_ν(T) = coef * Z² * g_ff / (sqrt(T) * ν³).

See `gaunt_ff_vanHoof` for details about where our gaunt factor data comes from. For simplicity, we
enforce temperature and ν bounds constraints that don't include all of the available gaunt factor
data. In practice, this should never be a concern for stellar spectroscopy.
"""
function hydrogenic_ff_absorption(ν::Real, T::Real, Z::Integer, ni::Real, ne::Real)
    inv_T = 1.0/T
    Z2 = Z*Z

    hν_div_kT = (hplanck_eV/kboltz_eV) * ν * inv_T
    log_u = log10(hν_div_kT)
    log_γ2 = log10((Rydberg_eV/kboltz_eV) * Z2 * inv_T)
    gaunt_ff = gaunt_ff_vanHoof(log_u, log_γ2)

    F_ν = 3.6919e8*gaunt_ff*Z2*sqrt(inv_T)/(ν*ν*ν)

    ni*ne*F_ν*(1-exp(-hν_div_kT))
end
