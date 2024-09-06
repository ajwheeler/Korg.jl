module Peach1970

using Interpolations: linear_interpolation
using ...Korg: Species, @species_str

coeffs = Dict{Species,Any}()

coeffs[species"He II"] = let
    # this comes from table III of Peach 1970 for neutral Helium

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

    linear_interpolation((T_vals, σ_vals), table_vals; extrapolation_bc=0)
end

coeffs[species"C II"] = let
    # this comes from table III of Peach 1970 for neutral Carbon

    # σ denotes the energy of the photon (in units of RydbergH*Zeff², where Zeff is the net charge
    # of the species participating in the interaction
    σ_vals = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]

    # The temperature (in K)
    T_vals = [4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0, 11000.0, 12000.0, 13000.0,
        14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 21000.0, 22000.0,
        23000.0, 24000.0, 25000.0, 26000.0, 27000.0, 28000.0, 29000.0, 30000.0, 32000.0,
        34000.0, 36000.0]

    # the (unitless) departure term
    table_vals = [-0.145 -0.144 -0.068 0.054 0.200 0.394;
                  -0.132 -0.124 -0.045 0.077 0.222 0.415;
                  -0.121 -0.109 -0.027 0.097 0.244 0.438;
                  -0.112 -0.095 -0.010 0.115 0.264 0.461;
                  -0.104 -0.082 0.005 0.133 0.284 0.484;
                  -0.095 -0.070 0.020 0.150 0.303 0.507;
                  -0.087 -0.058 0.034 0.166 0.321 0.529;
                  -0.079 -0.047 0.048 0.181 0.339 0.550;
                  -0.071 -0.036 0.061 0.196 0.356 0.570;
                  -0.063 -0.025 0.074 0.210 0.372 0.590;
                  -0.055 -0.015 0.086 0.223 0.388 0.609;
                  -0.047 -0.005 0.098 0.237 0.403 0.628;
                  -0.040 0.005 0.109 0.249 0.418 0.646;
                  -0.032 0.015 0.120 0.261 0.432 0.664;
                  -0.025 0.024 0.131 0.273 0.446 0.680;
                  -0.017 0.034 0.141 0.285 0.459 0.697;
                  -0.010 0.043 0.152 0.296 0.472 0.713;
                  -0.003 0.051 0.161 0.307 0.485 0.728;
                  0.004 0.060 0.171 0.317 0.497 0.744;
                  0.011 0.069 0.181 0.327 0.509 0.758;
                  0.018 0.077 0.100 0.337 0.521 0.773;
                  0.025 0.085 0.109 0.347 0.532 0.787;
                  0.032 0.093 0.208 0.356 0.543 0.800;
                  0.039 0.101 0.216 0.365 0.554 0.814;
                  0.046 0.109 0.225 0.374 0.564 0.827;
                  0.052 0.117 0.233 0.383 0.574 0.839;
                  0.059 0.124 0.241 0.391 0.585 0.852;
                  0.072 0.139 0.257 0.408 0.604 0.876;
                  0.085 0.154 0.273 0.424 0.623 0.900;
                  0.097 0.168 0.288 0.439 0.641 0.923]

    linear_interpolation((T_vals, σ_vals), table_vals; extrapolation_bc=0)
end

#coeffs[species"C III"] = let 
#    # this comes from table III of Peach 1970 for singly ionized Carbon. Peach broke this
#    # information up into 2 sub-tables:
#    # 1. data for a parent term ¹S, corresponding to the 1s²2s² ground state
#    # 2. data for a parent term ³Pᵒ, corresponding to the 1s²2s2p excited state
#
#    # Rather than deal with this complication in Korg's ff code, we omit the contribution from 
#    # the (higher energy) ³Pᵒ parent term
#
#    # σ denotes the energy of the photon (in units of RydbergH*Zeff², where Zeff is the net charge
#    # of the species participating in the interaction.
#    σ_vals = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30] # This is the same for both subtables
#
#    # First, specify the info for the 'S parent term:
#    # The temperature (in K)
#    T_vals_S = [16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 21000.0, 22000.0, 23000.0, 24000.0,
#                25000.0, 26000.0, 27000.0, 28000.0, 29000.0, 30000.0, 32000.0, 34000.0, 36000.0,
#                38000.0, 40000.0, 42000.0, 44000.0, 46000.0, 48000.0]
#
#    # the (unitless) departure term
#    table_vals_S = [0.032 0.135 0.247 0.346 0.432 0.516;
#                    0.034 0.136 0.246 0.344 0.429 0.511;
#                    0.036 0.137 0.246 0.343 0.426 0.507;
#                    0.038 0.138 0.246 0.342 0.424 0.504;
#                    0.040 0.140 0.247 0.341 0.423 0.502;
#                    0.042 0.141 0.247 0.341 0.422 0.500;
#                    0.044 0.142 0.248 0.341 0.421 0.499;
#                    0.046 0.144 0.249 0.341 0.420 0.497;
#                    0.048 0.145 0.249 0.341 0.420 0.496;
#                    0.050 0.146 0.250 0.341 0.419 0.495;
#                    0.051 0.148 0.251 0.341 0.419 0.495;
#                    0.053 0.149 0.252 0.342 0.419 0.494;
#                    0.055 0.150 0.252 0.342 0.418 0.493;
#                    0.056 0.151 0.253 0.342 0.418 0.492;
#                    0.058 0.153 0.254 0.342 0.418 0.492;
#                    0.061 0.155 0.255 0.343 0.418 0.491;
#                    0.063 0.157 0.257 0.343 0.417 0.490;
#                    0.066 0.159 0.258 0.343 0.417 0.488;
#                    0.068 0.161 0.259 0.344 0.416 0.487;
#                    0.071 0.163 0.260 0.344 0.416 0.486;
#                    0.073 0.165 0.261 0.344 0.415 0.485;
#                    0.075 0.167 0.262 0.344 0.414 0.483;
#                    0.078 0.168 0.262 0.344 0.413 0.482;
#                    0.080 0.169 0.263 0.343 0.412 0.481]
#
#    # First, specify the info for the 'S parent term:
#    # The temperature (in K)
#    # T_vals_P = [16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 21000.0, 22000.0, 23000.0, 24000.0,
#    #            25000.0, 26000.0, 27000.0, 28000.0, 29000.0, 30000.0, 32000.0, 34000.0, 36000.0,
#    #            38000.0, 40000.0, 42000.0, 44000.0, 46000.0, 48000.0]
#
#    # the (unitless) departure term
#    # table_vals_P = [0.031 0.161 0.319 0.478 0.631 0.809;
#    #                 0.033 0.162 0.318 0.475 0.627 0.803;
#    #                 0.035 0.164 0.318 0.474 0.624 0.799;
#    #                 0.038 0.165 0.318 0.473 0.622 0.797;
#    #                 0.040 0.166 0.319 0.472 0.621 0.795;
#    #                 0.042 0.168 0.319 0.472 0.621 0.794;
#    #                 0.044 0.169 0.320 0.472 0.620 0.793;
#    #                 0.046 0.171 0.321 0.473 0.620 0.793;
#    #                 0.048 0.173 0.322 0.473 0.621 0.793;
#    #                 0.050 0.174 0.323 0.474 0.621 0.793;
#    #                 0.051 0.176 0.324 0.475 0.622 0.793;
#    #                 0.053 0.177 0.325 0.475 0.622 0.794;
#    #                 0.055 0.179 0.327 0.476 0.623 0.795;
#    #                 0.057 0.180 0.328 0.477 0.624 0.796;
#    #                 0.059 0.182 0.329 0.478 0.625 0.797;
#    #                 0.062 0.185 0.332 0.480 0.627 0.799;
#    #                 0.066 0.188 0.334 0.483 0.629 0.801;
#    #                 0.069 0.191 0.337 0.485 0.631 0.804;
#    #                 0.072 0.194 0.339 0.487 0.633 0.806;
#    #                 0.075 0.197 0.342 0.489 0.635 0.809;
#    #                 0.078 0.200 0.345 0.491 0.638 0.811;
#    #                 0.081 0.203 0.347 0.493 0.640 0.814;
#    #                 0.084 0.206 0.349 0.406 0.642 0.816;
#    #                 0.087 0.209 0.352 0.498 0.644 0.819]
#
#    linear_interpolation((T_vals_S, σ_vals), table_vals_S, extrapolation_bc=0)
#end

coeffs[species"Si II"] = let
    # this comes from table III of Peach 1970 for neutral Silicon

    # σ denotes the energy of the photon (in units of RydbergH*Zeff², where Zeff is the net charge
    # of the species participating in the interaction
    σ_vals = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]

    # The temperature (in K)
    T_vals = [4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0, 11000.0, 12000.0, 13000.0,
        14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 21000.0, 22000.0,
        23000.0, 24000.0, 25000.0, 26000.0, 27000.0, 28000.0, 29000.0, 30000.0, 32000.0,
        34000.0, 36000.0]

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

    linear_interpolation((T_vals, σ_vals), table_vals; extrapolation_bc=0)
end

coeffs[species"Mg II"] = let
    # this comes from table III of Peach 1970 for neutral Magnesium

    # σ denotes the energy of the photon (in units of RydbergH*Zeff², where Zeff is the net charge
    # of the species participating in the interaction
    σ_vals = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]

    # The temperature (in K)
    T_vals = [4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000,
        17000, 18000, 19000, 20000, 21000, 22000, 23000, 24000, 25000, 26000, 27000, 28000,
        29000, 30000, 32000, 34000]

    table_vals = [-0.070 0.008 0.121 0.221 0.274 0.356;
                  -0.067 0.003 0.104 0.105 0.244 0.325;
                  -0.066 -0.002 0.091 0.175 0.221 0.302;
                  -0.065 -0.007 0.080 0.157 0.201 0.282;
                  -0.065 -0.012 0.069 0.141 0.183 0.264;
                  -0.065 -0.016 0.059 0.126 0.166 0.248;
                  -0.065 -0.020 0.049 0.113 0.151 0.232;
                  -0.066 -0.024 0.040 0.100 0.137 0.218;
                  -0.066 -0.028 0.032 0.088 0.124 0.205;
                  -0.066 -0.032 0.025 0.077 0.112 0.194;
                  -0.066 -0.035 0.018 0.067 0.101 0.183;
                  -0.066 -0.037 0.012 0.058 0.091 0.173;
                  -0.066 -0.040 0.006 0.049 0.082 0.164;
                  -0.066 -0.042 0.001 0.042 0.074 0.157;
                  -0.066 -0.044 -0.004 0.036 0.067 0.150;
                  -0.065 -0.045 -0.007 0.030 0.061 0.144;
                  -0.064 -0.046 -0.011 0.025 0.056 0.139;
                  -0.063 -0.047 -0.014 0.020 0.051 0.135;
                  -0.062 -0.048 -0.016 0.017 0.048 0.131;
                  -0.061 -0.048 -0.018 0.014 0.045 0.128;
                  -0.059 -0.047 -0.019 0.011 0.042 0.126;
                  -0.057 -0.047 -0.020 0.009 0.040 0.124;
                  -0.055 -0.046 -0.020 0.008 0.039 0.123;
                  -0.053 -0.045 -0.021 0.007 0.038 0.123;
                  -0.051 -0.044 -0.020 0.006 0.038 0.123;
                  -0.048 -0.042 -0.020 0.006 0.038 0.123;
                  -0.045 -0.040 -0.019 0.006 0.039 0.124;
                  -0.039 -0.035 -0.016 0.008 0.042 0.128;
                  -0.032 -0.030 -0.012 0.011 0.046 0.133]

    linear_interpolation((T_vals, σ_vals), table_vals; extrapolation_bc=0)
end

"""
    Peach1970.departure_coefficients()

This module contains interpolators of the tabulated ff departure coeffients from
[Peach+ 1970](https://ui.adsabs.harvard.edu/abs/1970MmRAS..73....1P/abstract), which we use to
correct the hydrogenic ff absorption coefficient for H I ff, C I ff, Si I ff, and Mg I ff.
It contains a dictionary (returned by `departure_coefficients()`), which maps `Species` to
interpolator objects.  Crucially, the dictionary is indexed by the species which actually
participates in the interaction, not the one after which the interaction is named.

Outside the regime in which Peach 1970 provides data, the interpolators return 0, falling back to
the hydreogenic approximation.

The species for which we use corrections are the same species which get corrected in
MARCS/Turbospectrum (see Table 1 of
[Gustafsson+ 2008](https://ui.adsabs.harvard.edu/abs/2008A%26A...486..951G/abstract)).
The choices seem are largely motivated by which species have departure terms at normal
stellar atmosphere conditions and which species are most abundant in the sun. For C II ff,
we include only the contribution from the ¹S parent term, even though (in contrast to other
speices) information is available for the ³Pᵒ term as well.

The free-free absorption coefficient (including stimulated emission) is given by:

``\\alpha_{\rm ff} = \\alpha_{\rm hydrogenic, ff}(\\nu, T, n_i, n_e; Z) (1 + D(T, \\sigma))``,

where

  - ``\\alpha_{\rm hydrogenic, ff}(\\nu, T, n_i, n_e; Z)`` should include the correction for
    stimulated emission.
  - ``n_i`` is the number density fo the ion species that participates in the interation, not the
    species the interaction is named after.
  - ``n_e`` is the number density of free electrons.
  - ``D(T, \\sigma)`` is specified as the `departure` arg, and is expected to interpolate over
    the tabulated values specified in Table III of Peach (1970).
  - σ denotes the energy of the photon in units of RydbergH*Zeff²

It might not be immediately obvious how the above equation relates to the equations presented in
Peach (1970). Peach describes the calculation for ``k_\nu^F``, the free-free absorption
coefficient (uncorrected for stimulated emission) per particle of the species that the interaction
is named after. In other words, he computes:

``k_\nu^F = \\alpha_{\rm ff}/n_{i-1} \\left(1 - e^\\frac{-h\\nu}{k T}\\right)^{-1}``,

where ``n_{i-1}`` is the number density of the species that the interaction is named after.
``k_\nu^F`` can directly be computed, under LTE, from just ``\nu``, ``T``, and ``n_{i-1}`` (the
Saha Equation relates ``\\alpha_{\rm ff}``'s dependence on ``n_e`` and ``n_i`` to ``n_{i-1}`` and
``T``.  Gray (2005) follows a similar convention when describing free-free absorption.

!!! warning


The tabulated data in this module was taken from Peach 1970 using OCR software, and may
contain mis-read values, although they produce reasonable behavior and there are no obvious
problems.
"""
departure_coefficients() = coeffs

end #module
