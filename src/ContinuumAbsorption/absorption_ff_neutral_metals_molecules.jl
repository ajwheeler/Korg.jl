using StaticArrays: SA, SMatrix
using Interpolations: LinearInterpolation, Throw

# This file defines functions to compute the free-free linear absorption coefficients for
# interactions involving neutral metals or neutral molecules. Based on convention, these
# interactions are all named after negative ions.

# -----------------------------------------------------------------------------------------
# Define absorption functions by combining data from both
# - [John 1975, MNRAS, 170, 5](https://ui.adsabs.harvard.edu/abs/1975MNRAS.170....5J)
# - [John 1975, MNRAS, 172, 305](https://ui.adsabs.harvard.edu/abs/1975MNRAS.172..305J)
# This is the procedure recommended by the latter paper.
# -----------------------------------------------------------------------------------------

"""
    _combined_john_neg_ion_absorption(ν, T, ndens_neutral, ne, short_wavelength_interp,
                                      long_wavelength_interp)

uses absorption data from John 1975, MNRAS, 172, 305 for 0.5 μm ≤ λ ≤ 10 μm and data from
John 1975, MNRAS, 170, 5 for λ > 10 μm. Returns linear absorption coefficient (in cm⁻¹)

These values are extremely uncertain.
"""
function _combined_john_neg_ion_absorption(ν, T, ndens_neutral, ne, short_wavelength_interp,
                                           long_wavelength_interp)
    # this involves a small amount of guessing, but we will confirm this is correct by comparing
    # He⁻ absorption coefficient with our reference

    λ_cm = c_cgs / ν # in cm
    Pₑ = ne * kboltz_cgs * T # partial pressure contributed by electrons

    cross_section = if λ_cm > 1e-3
        (λ_cm * λ_cm * 1e16) * long_wavelength_interp(T)
    else
        short_wavelength_interp(λ_cm, T)
    end

    ndens_neutral * Pₑ * cross_section
end


# load free-free absorption data from Table I of for interactions named after negative ions (but
# involve neutral species) John T. L. 1975, MNRAS, 170, 5. This data is for large λ (>~ 10 μm).
#
# MARCS/Turbospectrum only uses this data to compute for a small subset of species
const _long_wavelength_ff_interps = let
    T_vals = SA[ 50.0,    100.0,    500.0,    1000.0,   5000.0,   10000.0,   25000.0]

    table = Dict(
        "He"  => [    1.34,     0.968,    0.475,     0.353,    0.173,     0.13,      0.07  ],
        "Ne"  => [    0.107,    0.0823,   0.0604,    0.056,    0.041,     0.0328],  #-----
        "Ar"  => [    1.4,      0.673,    0.057,     0.019,    0.038,     0.062,     0.097 ],
        "Kr"  => [    6.17,     3.37,     0.527,     0.17,     0.034,     0.073,     0.138 ],
        "Xe"  => [   27.7,     13.8,      1.74,      0.548,    0.127,     0.234,     0.303 ],
        "Li"  => [   20.5,     41.6,     60.0,      28.0,      0.93,      0.062,     0.012 ],
        "Na"  => [   18.6,     24.6,     47.0,      26.0,      2.7,       1.3],     #-----
        "Cs"  => [ 1695.0,    853.0,     88.0,      57.0,      4.7,       2.8,       1.2   ],
        "Hg"  => [    5.09,     3.02,     0.845,     1.54],   #-----      -----      -----
        "N"   => [    0.263,    0.195,    0.13,      0.12,     0.1,       0.0904,    0.0699],
        "O"   => [    0.280,    0.202,    0.119,     0.103,    0.0743,    0.0654,    0.0512],
        "H2"  => [    1.9,      1.47,     0.866,     0.703,    0.429,     0.29,      0.112 ],
        "N2"  => [    0.597,    0.554,    0.506,     0.459,    0.333,     0.293,     0.14  ],
        "O2"  => [    0.762,    0.54,     0.343,     0.261,    0.171,     0.114,     0.076 ],
        "CO"  => [   21.3,      1.14,     0.642,     0.605,    0.45,      0.28,      0.15  ],
        "CO2" => [   43.5,     21.7,      4.1,       1.75,     0.192,     0.16,      0.13  ],
        "H20" => [  996.0,    342.0,     28.0]      #-----     -----      -----      -----
    )

    Dict(k => LinearInterpolation(view(T_vals, 1:length(v)), v.*1e-34, extrapolation_bc = Throw())
         for (k,v) in table)
end

# Most John1975b tables share a common set of λ and T values.
const _typical_john75b_wavelengths = SA[1.0e-3, 5.0e-4, 2.5e-4, 1.5e-4, 1.0e-4, 7.5e-5, 5.0e-5] #cm
const _typical_john75b_T_vals = SA[100.0, 500.0, 1000.0, 2500.0, 5000.0, 7500.0, 10000.0, 12500.0,
                                   15000.0] #K

function _build_interpolator_john75b(λ_vals, T_vals, elements)
    LinearInterpolation((reverse(λ_vals), T_vals),
                        view(elements,size(elements)[1]:-1:1,:), extrapolation_bc=0)
end

# for testing only!  We have a better He- ff treatment.
const _short_wavelength_Heminus_ff_interp = let
    # this comes from Table I of John 1975, MNRAS, 172, 305 and it holds He⁻ ff absorption
    # cross sections (including the correction for stimulated emission)
    elements = SMatrix{7,9}(
        [106.0  44.1  34.1  23.5   17.3   14.2   12.2   10.9    9.78;
          33.1  10.8   8.16  5.78   4.30   3.54   3.04   2.72   2.44;
          11.0   2.94  2.00  1.39   1.06   0.876  0.756  0.679  0.0609;
           4.98  1.21  0.756 0.484  0.370  0.310  0.269  0.243  0.218;
           2.67  0.620 0.366 0.214  0.159  0.134  0.118  0.107  0.0961;
           1.72  0.389 0.224 0.122  0.0880 0.0740 0.0650 0.0596 0.0535;
           0.931 0.205 0.114 0.0579 0.0388 0.0319 0.0281 0.0263 0.0233] .* 1e-26
    )

    _build_interpolator_john75b(_typical_john75b_wavelengths, _typical_john75b_T_vals, elements)
end


const _short_wavelength_Ominus_ff_interp = let
    # In the future, we may want to use the O⁻ data from John & Williams (1975),
    # [John & Williams (1977)](https://ui.adsabs.harvard.edu/abs/1977JQSRT..17..169J/abstract), or
    # [Chung & Lin 1995](https://ui.adsabs.harvard.edu/abs/1995PhRvA..51.1221C/abstract).
    #
    # Unfortunately, Chung & Lin 1995 only provide data for T ≥ 5000 K

    # this comes from Table I of John 1975, MNRAS, 172, 305 and it holds O⁻ ff absorption
    # cross sections (including the correction for stimulated emission)

    elements = SMatrix{7,9}(
        [20.1  10.9    9.86   8.38   7.41   6.91   6.54   6.22   5.97;
          6.26  2.61   2.32   2.05   1.84   1.72   1.63   1.55   1.49;
          2.05  0.693  0.556  0.489  0.450  0.425  0.405  0.387  0.371;
          0.924 0.281  0.207  0.160  0.157  0.150  0.144  0.138  0.133;
          0.494 0.142  0.0992 0.0739 0.0673 0.0647 0.0626 0.0603 0.0583;
          0.318 0.0885 0.0600 0.0420 0.0396 0.0355 0.0344 0.0333 0.0323;
          0.171 0.0461 0.0302 0.0196 0.0160 0.0152 0.0147 0.0143 0.0139] .* 1e-26
    )

    _build_interpolator_john75b(_typical_john75b_wavelengths, _typical_john75b_T_vals, elements)
end

"""
    Ominus_ff(ν, T, nH2_I, ne; kwargs...)

Compute the O⁻ free-free linear absorption coefficient (in cm⁻¹).

The naming scheme for free-free absorption is counter-inutitive. This actually refers to the
reaction:  `photon + e⁻ + O I -> e⁻ + O I.`

# Arguments
- `ν::AbstractVector{<:Real}`: sorted frequency vector in Hz
- `T`: temperature in K
- `nO_I`: the total number density of neutral O.
- `ne`: the number density of free electrons.

For a description of the kwargs, see [Continuum Absorption Kwargs](@ref).

# Notes

This function is defined for 100 K ≤ T ≤ 15000 K and 0.5 μ ≤ λ. It uses data from 
[John 1975, MNRAS, 170, 5](https://ui.adsabs.harvard.edu/abs/1975MNRAS.170....5J) for λ > 10 μm
and data from [John 1975, MNRAS, 172, 305](https://ui.adsabs.harvard.edu/abs/1975MNRAS.172..305J)
for λ ≤ 10 μm. This is the approach recomended by the latter paper. The error on these absorption
coefficients is probably large.
"""
Ominus_ff = bounds_checked_absorption(
    (ν, T, nO_I, ne) -> _combined_john_neg_ion_absorption(ν, T, nO_I, ne,
                                                           _short_wavelength_Ominus_ff_interp,
                                                           _long_wavelength_ff_interps["O"]),
    ν_bound = λ_to_ν_bound(closed_interval(5.0e-5, Inf)),
    temp_bound = closed_interval(100.0, 15000.0)
)

const _short_wavelength_H2minus_ff_interp = let
    # this comes from Table II of John 1975, MNRAS, 172, 305 and it holds H₂⁻ ff absorption
    # cross sections (including the correction for stimulated emission)
    #
    # I manually entered this table (so it's plausible that there are transcription errors)

    elements = SMatrix{7,9}(
        [153.0  79.7    67.9   54.6   42.7    34.8    29.3    24.6    20.9;
          47.4  19.3    16.2   13.4   10.6     8.69    7.29    6.13    5.21;
          15.6   5.22    3.92   3.19   2.59    2.14    1.80    1.52    1.30;
           7.02  2.13    1.46   1.09   0.894   0.750   0.638   0.542   0.463;
           3.76  1.07    0.695  0.470  0.380   0.322   0.277   0.237   0.203;
           2.41  0.688   0.419  0.265  0.207   0.176   0.152   0.130   0.112;
           1.30  0.347   0.210  0.123  0.0893  0.0748  0.0647  0.0558  0.0482] .* 1e-26
    )

    _build_interpolator_john75b(_typical_john75b_wavelengths, _typical_john75b_T_vals, elements)
end

"""
    H2minus_ff(ν, T, nH2_I, ne; kwargs...)

Compute the H₂⁻ free-free linear absorption coefficient (in cm⁻¹).

The naming scheme for free-free absorption is counter-inutitive. This actually refers to the
reaction:  `photon + e⁻ + H₂ I -> e⁻ + H₂ I.`

# Arguments
- `ν::AbstractVector{<:Real}`: sorted frequency vector in Hz
- `T`: temperature in K
- `nH2_I`: the total number density of neutral H₂.
- `ne`: the number density of free electrons.

For a description of the kwargs, see [Continuum Absorption Kwargs](@ref).

# Notes

This function is defined for 100 K ≤ T ≤ 15000 K and 0.5 μ ≤ λ. It uses data from 
[John 1975, MNRAS, 170, 5](https://ui.adsabs.harvard.edu/abs/1975MNRAS.170....5J) for λ > 10 μm
and data from [John 1975, MNRAS, 172, 305](https://ui.adsabs.harvard.edu/abs/1975MNRAS.172..305J)
for λ ≤ 10 μm. This is the approach recomended by the latter paper. The error on these absorption
coefficients is probably large.
"""
H2minus_ff = bounds_checked_absorption(
    (ν, T, nH2_I, ne) -> _combined_john_neg_ion_absorption(ν, T, nH2_I, ne,
                                                           _short_wavelength_H2minus_ff_interp,
                                                           _long_wavelength_ff_interps["H2"]),
    ν_bound = λ_to_ν_bound(closed_interval(5.0e-5, Inf)),
    temp_bound = closed_interval(100.0, 15000.0)
)


const _short_wavelength_COminus_ff_interp = let
    # this comes from Table II of John, T. L. 1975, MNRAS, 172, 305 and it holds CO⁻ ff absorption
    # cross sections (including the correction for stimulated emission) at smaller λ

    # the fact that the last column has numbers that are larger than preceding column isn't a typo
    elements = SMatrix{7,9}(
        [212.0  61.5  58.2   53.7   44.9    42.9    24.4    23.1    25.0;
          69.7  15.1  13.7   13.6   12.2    10.0     6.24    6.18    6.21;
          23.7  4.01   3.24   3.20   2.88    2.41    1.73    1.54    1.48;
          10.8  1.62   1.22   1.15   0.917   0.907   0.590   0.500   0.512;
           5.86 0.889  0.643  0.478  0.378   0.367   0.236   0.215   0.232;
           3.78 0.524  0.345  0.236  0.192   0.196   0.129   0.120   0.132;
           2.03 0.255  0.161  0.106  0.0830  0.0850  0.0557  0.0525  0.0596] .* 1e-26
    )

    _build_interpolator_john75b(_typical_john75b_wavelengths, _typical_john75b_T_vals, elements)
end

"""
    COminus_ff(ν, T, nCO_I, ne; kwargs...)

Compute the CO⁻ free-free linear absorption coefficient (in cm⁻¹).

The naming scheme for free-free absorption is counter-inutitive. This actually refers to the
reaction:  `photon + e⁻ + CO I -> e⁻ + CO I.`

# Arguments
- `ν::AbstractVector{<:Real}`: sorted frequency vector in Hz
- `T`: temperature in K
- `nCO_I`: the total number density of neutral CO.
- `ne`: the number density of free electrons.

For a description of the kwargs, see [Continuum Absorption Kwargs](@ref).

# Notes

This function is defined for 100 K ≤ T ≤ 15000 K and 0.5 μ ≤ λ. It uses data from 
[John 1975, MNRAS, 170, 5](https://ui.adsabs.harvard.edu/abs/1975MNRAS.170....5J) for λ > 10 μm
and data from [John 1975, MNRAS, 172, 305](https://ui.adsabs.harvard.edu/abs/1975MNRAS.172..305J)
for λ ≤ 10 μm. This is the approach recomended by the latter paper. The error on these absorption
coefficients is probably large.
"""
COminus_ff = bounds_checked_absorption(
    (ν, T, nCO_I, ne) -> _combined_john_neg_ion_absorption(ν, T, nCO_I, ne,
                                                           _short_wavelength_COminus_ff_interp,
                                                           _long_wavelength_ff_interps["CO"]),
    ν_bound = λ_to_ν_bound(closed_interval(5.0e-5, Inf)),
    temp_bound = closed_interval(100.0, 15000.0)
)


# set aside the following for the future
#const _small_wavelength_H2Ominus_ff_cross_section = let
#    # this comes from Table II of John 1975, MNRAS, 172, 305 and it holds H₂O⁻ ff absorption
#    # cross sections (including the correction for stimulated emission) at smaller λ
#
#    λ_vals = SA[10.0, 5.0, 2.5, 1.5, 1.0] .* 1e-6 # convert to cm
#    T_vals = SA[100.0, 500.0, 1000.0, 2500.0, 5000.0]
#
#    elements = SMatrix{5,5}([59160.0 4906.0 1247.0 163.0 25.0;
#                             19867.0 1539.0  394.0  53.0  8.2;
#                              6850.0  503.0  126.0  18.5  3.2;
#                              3153.0  227.0   55.2   8.4  NaN;
#                              1708.0  122.0   29.2   NaN  NaN] .* 1e-26)
#
#    _build_interpolator_john75b(λ_vals, T_vals, elements)
#end


# -----------------------------------------------------------------------------------------
# Define absorption functions using data from papers where Bell & Berrington were coauthors
# -----------------------------------------------------------------------------------------

"""
    _bell_berrington_negative_ion_ff_absorption(ν, T, ndens_neutral, ne, cross_section_interp)
    
Each paper presenting absorption coefficient data that Bell & Berrington are co-authors on all 
seem to present their data in a similar format
"""
function _bell_berrington_ff_neg_ion_absorption(ν, T, ndens_neutral, ne, cross_section_interp)
    # cross_section_interp should include the correction for stimulated emission
    θ = 5040.0 / T
    λ = c_cgs * 1.0e8 / ν # in ångstroms
    Pₑ = ne * kboltz_cgs * T # partial pressure contributed by electrons
    ndens_neutral * Pₑ * cross_section_interp(λ, θ)
end

const _Cminus_ff_coef = let
    # this data was transcribed from table 3 of Bell+ 1988
    # (https://ui.adsabs.harvard.edu/abs/1988JPhB...21.2319B)
    #
    # these cross-sections include the stimulated emission correction

    λ_vals = SA[151890.0, 113918, 91134, 45567, 30378, 22784, 18227, 15189, 13019, 11392, 10126,
                9113, 7595, 6510, 5696, 5063, 4557, 3645, 3038]
    θ_vals = SA[0.5, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.8]

    entries =   [ 3.95e2 4.62e2 6.13e2 7.83e2 9.68e2 1.35e3 1.73e3 2.37e3;
                  2.23e2 2.61e2 3.47e2 4.43e2 5.48e2 7.67e2 9.78e2 1.34e3;
                  1.43e2 1.68e2 2.23e2 2.86e2 3.53e2 4.94e2 6.29e2 8.57e2;
                  3.72e1 4.40e1 5.93e1 7.66e1 9.51e1 1.33e2 1.69e2 2.31e2;
                  1.74e1 2.08e1 2.83e1 3.67e1 4.54e1 6.22e1 7.71e1 1.00e2;
                  1.04e1 1.26e1 1.75e1 2.28e1 2.81e1 3.79e1 4.61e1 5.74e1;
                  7.20   8.80   1.24e1 1.61e1 1.99e1 2.66e1 3.19e1 3.88e1;
                  5.41   6.69   9.50   1.25e1 1.54e1 2.06e1 2.47e1 3.04e1;
                  4.32   5.39   7.73   1.02e1 1.25e1 1.68e1 2.02e1 2.52e1;
                  3.60   4.52   6.54   8.61   1.06e1 1.42e1 1.71e1 2.15e1;
                  3.12   3.96   5.80   7.69   9.56   1.27e1 1.52e1 1.88e1;
                  2.79   3.57   5.28   7.05   8.80   1.15e1 1.37e1 1.70e1;
                  2.29   2.97   4.40   5.87   7.27   9.77   1.18e1 1.45e1;
                  1.98   2.57   3.85   5.15   6.42   8.70   1.06e1 1.31e1;
                  1.75   2.29   3.45   4.63   5.79   7.88   9.62   1.20e1;
                  1.58   2.08   3.14   4.22   5.29   7.25   8.88   1.12e1;
                  1.45   1.91   2.90   3.92   4.92   6.77   8.31   1.05e1;
                  1.22   1.62   2.49   3.39   4.29   5.95   7.35   9.33;
                  1.07   1.43   2.21   3.04   3.86   5.38   6.67   8.50] .* 1e-27

    LinearInterpolation((reverse(λ_vals), θ_vals),
                        view(entries,size(entries)[1]:-1:1,:), extrapolation_bc=Throw())
end

"""
    Cminus_ff(ν, T, nC_I, ne; kwargs...)

Compute the C⁻ free-free linear absorption coefficient (in cm⁻¹).

The naming scheme for free-free absorption is counter-inutitive. This actually refers to the
reaction:  `photon + e⁻ + C I -> e⁻ + C I.`

# Arguments
- `ν::AbstractVector{<:Real}`: sorted frequency vector in Hz
- `T`: temperature in K
- `nC_I`: the total number density of C I.
- `ne`: the number density of free electrons.

For a description of the kwargs, see [Continuum Absorption Kwargs](@ref).

# Notes
This function is defined for 1800 K ≤ T ≤ 10080 K and 3038 Å ≤ λ ≤ 159890 Å. It uses data from
[Bell et al. 1988](https://ui.adsabs.harvard.edu/abs/1988JPhB...21.2319B).
"""
Cminus_ff = bounds_checked_absorption(
    (ν, T, nC_I, ne) -> _bell_berrington_ff_neg_ion_absorption(ν, T, nC_I, ne, _Cminus_ff_coef),
    ν_bound = λ_to_ν_bound(closed_interval(3.038e-5, 1.51890e-3)),
    temp_bound = closed_interval(1800.0, 10080.0)
)


const _Nminus_ff_coef = let
    # digitized table 3 from Ramsbottom et al. 1992
    # (https://iopscience.iop.org/article/10.1088/0953-4075/25/7/015
    #
    # these cross-sections include the stimulated emission correction
    #
    # OCR was misbehaving, but I think I got the transcription to be correct

    λ_vals = SA[159890.0, 113918, 91134, 45567, 30378, 22784, 18227, 15189, 13019, 11392, 10126,
                9113, 7595, 6510, 5696, 5063, 4557, 3645, 3038] # in ångstroms
    θ_vals = SA[0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.8]
    
    entries = SA[327.0  372.0  494.0  667.0 903.0 1220.0 1630.0 2140.0 2780.0 6860.0;
                 186.0  212.0  285.0  388.0 529.0  715.0  955.0 1250.0 1620.0 3900.0;
                 120.0  138.0  188.0  258.0 354.0  479.0  639.0  837.0 1080.0 2520.0;
                  33.4   40.2   59.6   87.8 126.0  175.0  236.0  309.0  394.0  854.0;
                  17.2   21.8   35.1   54.6  80.6  113.0  153.0  200.0  252.0  522.0;
                  11.6   15.3   26.4   42.3  63.4   89.6  121.0  156.0  197.0  396.0;
                   8.88  12.2   21.9   35.8  53.9   76.1  102.0  132.0  165.0  327.0;
                   7.37  10.4   19.2   31.6  47.5   66.9   89.5  115.0  144.0  283.0;
                   6.41   9.22  17.2   28.4  42.7   59.9   79.9  103.0  128.0  250.0;
                   5.78   8.44  15.9   26.2  39.3   55.0   73.2   93.7  116.0  227.0;
                   5.29   7.78  14.7   24.3  36.2   50.6   67.2   86.0  107.0  209.0;
                   4.89   7.23  13.7   22.5  33.5   46.7   62.0   79.3   98.5  192.0;
                   4.28   6.36  12.0   19.7  29.2   40.7   53.9   68.9   85.7  168.0;
                   3.92   5.84  11.0   17.9  27.0   37.0   49.0   62.6   77.8  152.0;
                   3.74   5.59  10.5   17.2  25.5   35.6   47.2   60.5   75.3  149.0;
                   3.46   5.2    9.8   16.0  24.3   33.3   43.2   54.7   71.0  133.0;
                   3.28   4.93   9.31  15.1  22.5   31.8   41.6   52.5   67.8  127.0;
                   2.99   4.41   8.48  13.7  20.8   28.8   38.0   47.5   60.9  116.0;
                   2.89   4.32   8.42  13.7  20.4   28.7   37.5   47.6   60.3  113.0] .* 1e-27

    LinearInterpolation((reverse(λ_vals), θ_vals),
                        view(entries,size(entries)[1]:-1:1,:), extrapolation_bc=Throw())
end


"""
    Nminus_ff(ν, T, nN_I, ne; kwargs...)

Compute the N⁻ free-free linear absorption coefficient (in cm⁻¹).

The naming scheme for free-free absorption is counter-inutitive. This actually refers to the
reaction:  `photon + e⁻ + N I -> e⁻ + N I.`

# Arguments
- `ν::AbstractVector{<:Real}`: sorted frequency vector in Hz
- `T`: temperature in K
- `nN_I`: the total number density of N I.
- `ne`: the number density of free electrons.

For a description of the kwargs, see [Continuum Absorption Kwargs](@ref).

# Notes

This function is defined for 1800 K ≤ T ≤ 10080 K and 3038 Å ≤ λ ≤ 159890 Å. It uses data from
[Ramsbottom et al. 1992](https://iopscience.iop.org/article/10.1088/0953-4075/25/7/015).
"""
Nminus_ff = bounds_checked_absorption(
    (ν, T, nN_I, ne) -> _bell_berrington_ff_neg_ion_absorption(ν, T, nN_I, ne, _Nminus_ff_coef),
    ν_bound = λ_to_ν_bound(closed_interval(3.038e-5, 1.51890e-3)),
    temp_bound = closed_interval(1800.0, 10080.0)
)
