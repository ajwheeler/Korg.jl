"""
According to Gray (2005), the bound-free contributions from He⁻ are usually assumed to be 
negligible because it only has one bound level with an ionization energy 19 eV. Supposedly the 
population of that level is too small to be worth considering.

We are currently missing free-free and bound free contributions from He I.
"""
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
        throw(DomainError(n, " Unknown excited state properties for He I"))
    end
    nsdens_div_partition * g_n * exp(-energy_level / (kboltz_eV * T))
end

const _Heminus_ff_absorption_interp = let
    #OCR'd from John (1994) https://ui.adsabs.harvard.edu/abs/1994MNRAS.269..871J
    theta_ff_absorption_interp = [0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.8, 3.6]
    lambda_ff_absorption_interp = 1e4 .* [0.5063, 0.5695, 0.6509, 0.7594, 0.9113, 1.1391, 1.5188,
        1.8225, 2.2782, 3.0376, 3.6451, 4.5564, 6.0751, 9.1127, 11.390, 15.1878]

    ff_absorption = [0.033 0.036 0.043 0.049 0.055 0.061 0.066 0.072 0.078 0.100 0.121
                     0.041 0.045 0.053 0.061 0.067 0.074 0.081 0.087 0.094 0.120 0.145
                     0.053 0.059 0.069 0.077 0.086 0.094 0.102 0.109 0.117 0.148 0.178
                     0.072 0.079 0.092 0.103 0.114 0.124 0.133 0.143 0.152 0.190 0.227
                     0.102 0.113 0.131 0.147 0.160 0.173 0.186 0.198 0.210 0.258 0.305
                     0.159 0.176 0.204 0.227 0.247 0.266 0.283 0.300 0.316 0.380 0.444
                     0.282 0.311 0.360 0.400 0.435 0.466 0.495 0.522 0.547 0.643 0.737
                     0.405 0.447 0.518 0.576 0.625 0.670 0.710 0.747 0.782 0.910 1.030
                     0.632 0.698 0.808 0.899 0.977 1.045 1.108 1.165 1.218 1.405 1.574
                     1.121 1.239 1.435 1.597 1.737 1.860 1.971 2.073 2.167 2.490 2.765
                     1.614 1.783 2.065 2.299 2.502 2.681 2.842 2.990 3.126 3.592 3.979
                     2.520 2.784 3.226 3.593 3.910 4.193 4.448 4.681 4.897 5.632 6.234
                     4.479 4.947 5.733 6.387 6.955 7.460 7.918 8.338 8.728 10.059 11.147
                     10.074 11.128 12.897 14.372 15.653 16.798 17.838 18.795 19.685 22.747 25.268
                     15.739 17.386 20.151 22.456 24.461 26.252 27.882 29.384 30.782 35.606 39.598
                     27.979 30.907 35.822 39.921 43.488 46.678 49.583 52.262 54.757 63.395 70.580]
    linear_interpolation((lambda_ff_absorption_interp, theta_ff_absorption_interp), ff_absorption;
                         extrapolation_bc=Throw())
end

function _Heminus_ff(ν::Real, T::Real, nHe_I_div_partition::Real, ne::Real)
    λ = c_cgs * 1.0e8 / ν # Å
    θ = 5040.0 / T

    # K includes contribution from stimulated emission
    K = 1e-26 * _Heminus_ff_absorption_interp(λ, θ) # [cm^4/dyn]

    Pe = ne * kboltz_cgs * T # partial pressure contributed by electrons
    nHe_I_gs = ndens_state_He_I(1, nHe_I_div_partition, T)

    K * nHe_I_gs * Pe
end

"""
    Heminus_ff(ν, T, nHe_I_div_partition, ne; kwargs...)

Compute the He⁻ free-free opacity κ.

The naming scheme for free-free absorption is counter-inutitive. This actually refers to the
reaction:  `photon + e⁻ + He I -> e⁻ + He I.`

# Arguments

  - `ν::AbstractVector{<:Real}`: sorted frequency vector in Hz
  - `T`: temperature in K
  - `nHe_I_div_partition`: the total number density of H I divided by its partition function.
  - `ne`: the number density of free electrons.

For a description of the kwargs, see [Continuum Absorption Kwargs](@ref).

# Notes

This uses the tabulated values from
[John (1994)](https://ui.adsabs.harvard.edu/abs/1994MNRAS.269..871J/abstract).  The quantity K is
the same used by [Bell and Berrington (1987)](https://doi.org/10.1088/0022-3700/20/4/019).  See
[`Hminus_ff`](@ref) for an explanation.

According to John (1994), improved calculations are unlikely to alter the tabulated data for
λ > 1e4Å, "by more than about 2%." The errors introduced by the approximations for
5.06e3 Å ≤ λ ≤ 1e4 Å "are expected to be well below 10%."
"""
Heminus_ff = bounds_checked_absorption(_Heminus_ff;
                                       ν_bound=λ_to_ν_bound(closed_interval(5.063e-5, 1.518780e-03)),
                                       temp_bound=closed_interval(1400, 10080))
