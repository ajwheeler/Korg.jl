#=
    absorption_mol_photodissociation.jl

    CH and OH molecular photodissociation continuum opacity.

    Ported from SYNTHE/ATLAS12 CHOP and OHOP functions (Kurucz, priv. comm.;
    implemented by Castelli 2005a). Cross-section × partition function tables
    from crossch.dat and crossoh.dat.

    The tables store log10(σ·Q) on a 2D grid of (temperature, photon energy),
    where Q is the molecular partition function. A separate partition function
    Q(T) is interpolated independently and divided out to recover σ(T, E).
    The final opacity contribution is:

        α_CH = n(CH) · σ_CH(T, E) · (1 - exp(-hν/kT))

    and similarly for OH.

    The energy grid is:
        CH: E = 0.1, 0.2, ..., 10.5 eV  (105 bins, but only bins 20–105 used, i.e. E ≥ 2.0 eV)
        OH: E = 2.1, 2.2, ..., 15.0 eV  (130 bins)

    The temperature grid for cross-sections is:
        T = 2000, 2500, ..., 9000 K  (15 points, 500 K steps)

    The partition function grid is:
        T = 1000, 1200, ..., 9000 K  (41 points, 200 K steps)

    References:
        Kurucz, R.L. (private communication)
        Castelli, F. 2005a, Memorie della Società Astronomica Italiana Supplementi, 8, 44
=#

using DelimitedFiles

# ============================================================================
# CH photodissociation
# ============================================================================

# CH partition function: T = 1000 to 9000 K in 200 K steps (41 points)
# From SYNTHE CHOP routine
const _CH_partition_function_table = [
    203.741,  249.643,  299.341,  353.477,  412.607,  477.237,  547.817,
    624.786,  708.543,  799.463,  897.912, 1004.227, 1118.738, 1241.761,
   1373.588, 1514.481, 1664.677, 1824.394, 1993.801, 2173.050, 2362.234,
   2561.424, 2770.674, 2989.930, 3219.204, 3458.378, 3707.355, 3966.005,
   4234.155, 4511.604, 4798.135, 5093.554, 5397.593, 5709.948, 6030.401,
   6358.646, 6694.379, 7037.313, 7387.147, 7743.579, 8106.313
]

# Load the CH cross-section × partition function table from crossch.dat
# Returns a 15×105 matrix of log10(σ·Q) values
function _load_crossch(datadir::String)
    filepath = joinpath(datadir, "crossch.dat")
    # skip the 5 header lines, read 105 rows × 15 columns
    data = readdlm(filepath, skipstart=5)
    # data is 105×15 (rows=energy bins, cols=temperatures)
    # transpose to match the Fortran convention: CROSSCH(IT, N)
    return Matrix{Float64}(data')  # now 15×105
end

const crossch = _load_crossch(_data_dir)

"""
    _interp_CH_partition(T)

Interpolate the CH partition function at temperature T (K).
Linear interpolation on a 200 K grid from 1000 to 9000 K.
"""
function _interp_CH_partition(T::Real)
    if T >= 9000.0
        return _CH_partition_function_table[end]
    end
    # index into the 200 K grid (T = 1000, 1200, ..., 9000 K)
    it = clamp(floor(Int, (T - 1000.0) / 200.0) + 1, 1, 40)
    T_node = 800.0 + it * 200.0
    frac = (T - T_node) / 200.0
    return _CH_partition_function_table[it] + 
           (_CH_partition_function_table[it+1] - _CH_partition_function_table[it]) * frac
end

"""
    _CH_cross_section(crossch_table, T, E_eV)

Compute the CH photodissociation cross-section σ (cm²) at temperature T (K)
and photon energy E_eV (eV) by bilinear interpolation of the cross-section ×
partition function table, then dividing out the partition function.

Returns 0.0 if T ≥ 9000 K or E_eV is outside the valid range (2.0–10.5 eV).
"""
function _CH_cross_section(crossch_table::Matrix{Float64}, T::Real, E_eV::Real)
    if T >= 9000.0
        return zero(promote_type(typeof(T), typeof(E_eV)))
    end

    # Energy index: E = 0.1, 0.2, ..., 10.5 eV → index = E/0.1
    # Only valid for indices 20–104 (E = 2.0–10.4 eV, interpolating to 10.5)
    n = floor(Int, E_eV * 10.0)
    if n < 20 || n >= 105
        return zero(promote_type(typeof(T), typeof(E_eV)))
    end
    E_node = n * 0.1
    frac_E = (E_eV - E_node) / 0.1

    # Temperature index for cross-section table: T = 2000, 2500, ..., 9000 K
    it = clamp(floor(Int, (T - 2000.0) / 500.0) + 1, 1, 14)
    T_node = 1500.0 + it * 500.0
    frac_T = (T - T_node) / 500.0

    # Bilinear interpolation of log10(σ·Q)
    log_σQ_low  = crossch_table[it, n]   + (crossch_table[it, n+1]   - crossch_table[it, n])   * frac_E
    log_σQ_high = crossch_table[it+1, n] + (crossch_table[it+1, n+1] - crossch_table[it+1, n]) * frac_E
    log_σQ = log_σQ_low + (log_σQ_high - log_σQ_low) * frac_T

    # Convert from log10(σ·Q) to σ by dividing out the partition function
    Q = _interp_CH_partition(T)
    σQ = exp(log_σQ * log(10))  # 10^(log10(σ·Q))
    return σQ / Q  # cross-section in cm², since the table already encodes cm² units
end


# ============================================================================
# OH photodissociation
# ============================================================================

# OH partition function: T = 1000 to 9000 K in 200 K steps (41 points)
# From SYNTHE OHOP routine
const _OH_partition_function_table = [
    145.979,  178.033,  211.618,  247.053,  284.584,  324.398,  366.639,
    411.425,  458.854,  509.012,  561.976,  617.823,  676.626,  738.448,
    803.363,  871.437,  942.735, 1017.330, 1095.284, 1176.654, 1261.510,
   1349.898, 1441.875, 1537.483, 1636.753, 1739.733, 1846.434, 1956.883,
   2071.080, 2189.029, 2310.724, 2436.155, 2565.283, 2698.103, 2834.571,
   2974.627, 3118.242, 3265.366, 3415.912, 3569.837, 3727.077
]

# Load the OH cross-section × partition function table from crossoh.dat
function _load_crossoh(datadir::String)
    filepath = joinpath(datadir, "crossoh.dat")
    data = readdlm(filepath, skipstart=5)
    return Matrix{Float64}(data')  # 15×130
end

const crossoh = _load_crossoh(_data_dir)

"""
    _interp_OH_partition(T)

Interpolate the OH partition function at temperature T (K).
Linear interpolation on a 200 K grid from 1000 to 9000 K.
"""
function _interp_OH_partition(T::Real)
    if T >= 9000.0
        return _OH_partition_function_table[end]
    end
    it = clamp(floor(Int, (T - 1000.0) / 200.0) + 1, 1, 40)
    T_node = 800.0 + it * 200.0
    frac = (T - T_node) / 200.0
    return _OH_partition_function_table[it] + 
           (_OH_partition_function_table[it+1] - _OH_partition_function_table[it]) * frac
end

"""
    _OH_cross_section(crossoh_table, T, E_eV)

Compute the OH photodissociation cross-section σ (cm²) at temperature T (K)
and photon energy E_eV (eV) by bilinear interpolation, then dividing out
the partition function.

Returns 0.0 if T ≥ 9000 K or E_eV is outside the valid range (2.1–15.0 eV).
"""
function _OH_cross_section(crossoh_table::Matrix{Float64}, T::Real, E_eV::Real)
    if T >= 9000.0
        return zero(promote_type(typeof(T), typeof(E_eV)))
    end

    # Energy index: offset by 2.0 eV. E = 2.1, 2.2, ..., 15.0 eV
    # n = (E - 2.0) / 0.1, valid for n = 1..129
    n = floor(Int, (E_eV - 2.0) / 0.1)
    if n <= 0 || n >= 130
        return zero(promote_type(typeof(T), typeof(E_eV)))
    end
    E_node = n * 0.1 + 2.0
    frac_E = (E_eV - E_node) / 0.1

    # Temperature index: T = 2000, 2500, ..., 9000 K
    it = clamp(floor(Int, (T - 2000.0) / 500.0) + 1, 1, 14)
    T_node = 1500.0 + it * 500.0
    frac_T = (T - T_node) / 500.0

    # Bilinear interpolation of log10(σ·Q)
    log_σQ_low  = crossoh_table[it, n]   + (crossoh_table[it, n+1]   - crossoh_table[it, n])   * frac_E
    log_σQ_high = crossoh_table[it+1, n] + (crossoh_table[it+1, n+1] - crossoh_table[it+1, n]) * frac_E
    log_σQ = log_σQ_low + (log_σQ_high - log_σQ_low) * frac_T

    Q = _interp_OH_partition(T)
    σQ = exp(log_σQ * log(10))
    return σQ / Q
end


# ============================================================================
# Public interface: add CH + OH photodissociation opacity to α
# ============================================================================

"""
    mol_photodissociation_absorption!(α, νs, T, number_densities, crossch, crossoh)

Add the CH and OH molecular photodissociation continuum opacity to the
absorption coefficient vector α (cm⁻¹).

The cross-section tables `crossch` and `crossoh` should be loaded once at
startup via `_load_crossch` and `_load_crossoh`.

This function requires that Korg's chemical equilibrium solver has computed
number densities for CH (`species"CH"`) and OH (`species"OH"`).

# Arguments
  - `α`: absorption coefficient vector to be modified in-place (cm⁻¹)
  - `νs`: frequency vector (Hz), sorted
  - `T`: temperature (K)
  - `number_densities`: Dict mapping Species → number density (cm⁻³)

# References
  - Kurucz, R.L. (private communication)
  - Castelli, F. 2005a, Mem. S.A.It. Suppl., 8, 44

# Notes
  The cross-section tables encode log10(σ·Q) where Q is the molecular
  partition function. The Fortran CHOP/OHOP functions in SYNTHE divide
  out Q from the table and then multiply by n(molecule)/U(molecule) from
  the population arrays. Here we use the total number density n(CH) or
  n(OH) from Korg's chemical equilibrium, which already accounts for the
  partition function. Therefore we divide σ·Q by Q(T) to get σ, then
  multiply by n(molecule).

  In SYNTHE, the final COOLOP assembly multiplies CHOP(J) by:
      XNFP(J, 846) × STIM(J) / RHO(J)
  where XNFP(J, 846) = n(CH)/U(CH) and STIM = (1 - exp(-hν/kT)).
  Here, Korg's number_densities already give n(CH) directly, so we
  multiply by n(CH) · σ(T,E) · (1 - exp(-hν/kT)).
"""
function mol_photodissociation_absorption!(α::AbstractVector, νs::AbstractVector,
                                            T::Real, number_densities::Dict)
    # Physical constants
    h_eV = 4.135667696e-15  # Planck constant in eV·s
    k_eV = 8.617333262e-5   # Boltzmann constant in eV/K

    # println("Species in number_densities: ", keys(number_densities))

    # CH contribution
    nCH = get(number_densities, species"CH", 0.0)
    # println("n(CH) = ", nCH)
    if nCH > 0.0
        for i in eachindex(νs)
            E_eV = h_eV * νs[i]
            σ = _CH_cross_section(crossch, T, E_eV)
            if σ > 0.0
                stim = 1.0 - exp(-E_eV / (k_eV * T))
                α[i] += nCH * σ * stim
            end
        end
    end

    # OH contribution
    nOH = get(number_densities, species"OH", 0.0)
    # println("n(OH) = ", nOH)
    if nOH > 0.0
        for i in eachindex(νs)
            E_eV = h_eV * νs[i]
            σ = _OH_cross_section(crossoh, T, E_eV)
            if σ > 0.0
                stim = 1.0 - exp(-E_eV / (k_eV * T))
                α[i] += nOH * σ * stim
            end
        end
    end

    return α
end