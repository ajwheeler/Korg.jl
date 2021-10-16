using Interpolations: LinearInterpolation, Throw

"""
    metal_ff_absorption_departure(ν, T, Z, ni, ne, departure)

Computes the free-free linear absorption coefficient (in cm⁻¹) of a metal species for which a 
departure term is specified.

A free-free interaction is named as though the species interacting with the free electron had one 
more bound electron (in other words it's named as though the free-electron and ion were bound 
together). The `Z` argument and `ni` values should be specified for the species that actually
participates in the reaction. For completeness 2 examples provided below:
- Si I ff absorption: `ni` holds the number density of Si II, and `Z=1` (net charge of Si II)
- Si II ff absorption: `ni` holds the number density of Si III, and `Z=2` (net charge of Si III)

# Arguments
- `ν`: frequency in Hz
- `T`: temperature in K
- `Z::Integer`: the net charge of the ion species (that participates in the interaction).
- `ni`: the number density of the ion species (that participates in the interaction) in cm⁻³.
- `ne`: the number density of free electrons.
- `departure`: a function or functor that gives the departure (i.e. relative difference from
  hydrogenic absorption) value at a given Temperature and σ, where 
  σ = (ν / Hz) / (Z² Rydberg_eV / hplanck_eV).

# Notes

This calculation is based on the approach described in
[Peach 1970](https://ui.adsabs.harvard.edu/abs/1970MmRAS..73....1P/abstract).
In short, the free-free absorption (including stimulated emission) is given by:

``\\alpha_{\rm ff} = \\alpha_{\rm hydrogenic, ff}(\\nu, T, n_i, n_e; Z) (1 + D(T, \\sigma))``,

where
- ``\\alpha_{\rm hydrogenic, ff}(\\nu, T, n_i, n_e; Z)`` should includes the correction for 
  stimulated emission.
- ``D(T, \\sigma)`` is specified as the `departure` arg, and is expected to interpolate over
the tabulated values specified in Table III of Peach (1970).

It might not be immediately obvious the above equation relates to the equations presented in
Peach (1970). Peach describes the calculation for ``k_\nu^F``, the free-free absorption 
coefficient (uncorrected for stimulated emission) per particle of the species that the interaction 
is named after. In other words, he computes:
 
``k_\nu^F = \\alpha_{\rm ff}/n_{i-1} \\left(1 - e^\\frac{-h\\nu}{k T}\\right)^{-1}``,

where ``n_{i-1}`` is the number density of the species that the interaction is named after.
``k_\nu^F`` can directly be computed, under LTE, from just ``\nu``, ``T``, and ``n_{i-1}`` (the
Saha Equation relates ``\\alpha_{\rm ff}``'s dependence on ``n_e`` and ``n_i`` to ``n_{i-1}`` and 
``T``.  Gray (2005) follows a similar convention when describing free-free absorption.
"""
function metal_ff_absorption_departure(ν::Real, T::Real, Z::Integer, ni::Real, ne::Real, departure)
    σ = ν/ Z^2 * (hplanck_eV / Rydberg_eV)
    hydrogenic_ff_absorption(ν, T, Z, ni, ne) * (1 + departure(T, σ))
end


# should we be loading this from a hdf5 table?
const _departure_term_He_I_ff = begin
    # this comes from table III of Peach 1970 for neutral Helium
    # Warning: it's plausible that the OCR software made transcription errors

    # σ denotes the energy of the photon (in units of RydbergH*Zeff², where Zeff is the net charge
    # of the species participating in the interaction
    σ_vals = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]

    # The temperature (in K)
    T_vals = [10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0,
              19000.0, 20000.0, 21000.0, 22000.0, 23000.0, 24000.0, 25000.0, 26000.0, 27000.0,
              28000.0, 29000.0, 30000.0, 32000.0, 34000.0, 36000.0, 38000.0, 40000.0, 42000.0,
              44000.0, 46000.0, 48000.0]

    # the (unitless) departure term
    table_vals = [0.016 0.039 0.069 0.100 0.135 0.169;
                  0.018 0.041 0.071 0.103 0.137 0.172;
                  0.020 0.043 0.073 0.105 0.139 0.174;
                  0.022 0.045 0.075 0.107 0.142 0.176;
                  0.024 0.047 0.078 0.109 0.144 0.179;
                  0.026 0.050 0.080 0.112 0.146 0.181;
                  0.028 0.052 0.082 0.114 0.148 0.183;
                  0.029 0.054 0.084 0.116 0.151 0.185;
                  0.031 0.056 0.086 0.118 0.153 0.187;
                  0.033 0.058 0.088 0.120 0.155 0.190;
                  0.035 0.060 0.090 0.122 0.157 0.192;
                  0.037 0.062 0.092 0.125 0.159 0.194;
                  0.039 0.064 0.095 0.127 0.162 0.196;
                  0.041 0.066 0.097 0.129 0.164 0.198;
                  0.043 0.068 0.099 0.131 0.166 0.201;
                  0.045 0.070 0.101 0.133 0.168 0.203;
                  0.047 0.072 0.103 0.135 0.170 0.205;
                  0.049 0.074 0.105 0.138 0.173 0.207;
                  0.050 0.076 0.107 0.140 0.175 0.209;
                  0.052 0.079 0.109 0.142 0.177 0.211;
                  0.054 0.081 0.111 0.144 0.179 0.214;
                  0.058 0.085 0.115 0.148 0.183 0.218;
                  0.062 0.089 0.119 0.153 0.188 0.222;
                  0.065 0.093 0.124 0.157 0.102 0.226;
                  0.069 0.096 0.128 0.161 0.196 0.230;
                  0.072 0.100 0.132 0.165 0.200 0.235;
                  0.076 0.104 0.135 0.169 0.204 0.239;
                  0.079 0.108 0.139 0.173 0.208 0.243;
                  0.082 0.111 0.143 0.177 0.212 0.247;
                  0.085 0.115 0.147 0.181 0.216 0.251]

    LinearInterpolation((T_vals, σ_vals), table_vals, extrapolation_bc=Throw())
end

_He_I_ff(ν::Real, T::Real, ndens_HII::Real, nₑ::Real) =
    metal_ff_absorption_departure(ν, T, 1, ndens_HII, nₑ, _departure_term_He_I_ff)


const _departure_term_Si_I_ff = begin
    # this comes from table III of Peach 1970 for neutral Silicon
    # Warning: it's plausible that the OCR software made transcription errors

    # σ denotes the energy of the photon (in units of RydbergH*Zeff², where Zeff is the net charge
    # of the species participating in the interaction
    σ_vals = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]

    # The temperature (in K)
    T_vals = [4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000,
              17000, 18000, 19000, 20000, 21000, 22000, 23000, 24000, 25000, 26000, 27000, 28000,
              29000, 30000, 32000, 34000, 36000]
    
    table_vals = [-0.079 0.033 0.214 0.434 0.650 0.973;
                  -0.066 0.042 0.216 0.429 0.642 0.062;
                  -0.056 0.050 0.220 0.430 0.643 0.965;
                  -0.048 0.057 0.224 0.433 0.648 0.974;
                  -0.040 0.063 0.229 0.436 0.653 0.081;
                  -0.033 0.069 0.233 0.440 0.659 0.995;
                  -0.027 0.074 0.238 0.444 0.666 1.007;
                  -0.021 0.080 0.242 0.448 0.672 1.019;
                  -0.015 0.085 0.246 0.452 0.679 1.031;
                  -0.010 0.089 0.250 0.456 0.685 1.042;
                  -0.004 0.094 0.254 0.459 0.692 1.054;
                  0.001 0.009 0.258 0.463 0.698 1.065;
                  0.006 0.103 0.262 0.467 0.705 1.076;
                  0.011 0.107 0.265 0.471 0.711 1.087;
                  0.016 0.112 0.269 0.474 0.717 1.097;
                  0.021 0.116 0.273 0.478 0.724 1.108;
                  0.026 0.120 0.277 0.482 0.730 1.118;
                  0.030 0.125 0.281 0.486 0.736 1.127;
                  0.035 0.129 0.285 0.490 0.742 1.137;
                  0.040 0.134 0.289 0.493 0.747 1.146;
                  0.045 0.138 0.293 0.497 0.753 1.155;
                  0.050 0.143 0.297 0.501 0.759 1.164;
                  0.055 0.147 0.301 0.505 0.765 1.173;
                  0.060 0.152 0.305 0.509 0.770 1.181;
                  0.065 0.156 0.310 0.513 0.776 1.189;
                  0.071 0.161 0.314 0.517 0.781 1.197;
                  0.076 0.166 0.318 0.520 0.787 1.205;
                  0.087 0.176 0.328 0.528 0.798 1.221;
                  0.008 0.186 0.317 0.537 0.809 1.236;
                  0.109 0.196 0.346 0.545 0.819 1.251]

    LinearInterpolation((T_vals, σ_vals), table_vals, extrapolation_bc=Throw())
end
