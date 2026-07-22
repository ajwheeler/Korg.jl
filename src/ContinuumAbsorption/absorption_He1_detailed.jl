#=
    absorption_He1_detailed.jl

    Detailed He I bound-free opacity, ported from SYNTHE/ATLAS12 HE1OP.

    This ADDS He I bound-free opacity to Korg, which previously had NONE.
    (The only He-related continuum in Korg was He⁻ ff from John 1994.)

    Components:
      1. Cross-section functions for 5 resolved low-lying states
      2. Hydrogenic cross-sections for n=3 shell (5 states) and n=4–27
      3. Boltzmann population calculations for each level
      4. Inner-shell ionization (two-electron transitions)
      5. Dissolved-level contribution near the series limit

    He I free-free (He II + e⁻ bremsstrahlung) is NOT included here because
    it is already handled by positive_ion_ff_absorption! in Korg's
    total_continuum_absorption (via hydrogenic_ff_absorption with Z=2).

    == Integration ==

    In ContinuumAbsorption.jl, add the following call to total_continuum_absorption:

        # He I detailed bound-free (HE1OP port from SYNTHE)
        He1_detailed_bf!(α, νs, T, number_densities, partition_funcs)

    This should be placed after the Heminus_ff call and before
    positive_ion_ff_absorption!.

    == References ==
        Kurucz, R.L. 1970, SAO Special Report 309
        Karzas, W.J. & Latter, R. 1961, ApJS 6, 167
        Mathisen, R. 1984, Inst. Theor. Astrophys. Oslo Pub. No. 1
        Fernley, J.A., Taylor, K.T. & Seaton, M.J. 1987, J.Phys.B 20, 6457
=#


# ============================================================================
# Physical constants (prefixed to avoid collisions with Korg's own constants)
# ============================================================================

const _He1_h_eV       = 4.135667696e-15     # Planck constant (eV·s)
const _He1_k_eV       = 8.617333262e-5      # Boltzmann constant (eV/K)
const _He1_h_cgs      = 6.62607015e-27      # Planck constant (erg·s)
const _He1_k_cgs      = 1.380649e-16        # Boltzmann constant (erg/K)
const _He1_c_cgs      = 2.99792458e10       # speed of light (cm/s)
const _He1_c_ang      = 2.99792458e18       # speed of light (Å/s)
const _He1_Rydberg_He = 438908.8637         # He I Rydberg constant (cm⁻¹)
const _He1_ionization_eV = 24.587           # He I ionization potential (eV)
const _He1_dissolved_eV  = 23.730           # dissolved level threshold (eV)
const _He1_dissolved_coeff = 2.815e29       # dissolved-level σ coefficient


# ============================================================================
# Level metadata
# ============================================================================

const _He1_level_χ_eV = [0.0, 19.819, 20.615, 20.964, 21.217,
                          22.718, 22.920, 23.006, 23.073, 23.086]

const _He1_level_ν_threshold = [5.945209e15, 1.152844e15, 0.9603331e15,
                                 0.8761076e15, 0.8147104e15, 0.4519048e15,
                                 0.4030971e15, 0.3821191e15, 0.3660215e15,
                                 0.3627891e15]

const _He1_level_g = [1.0, 3.0, 1.0, 9.0, 3.0, 3.0, 1.0, 9.0, 20.0, 3.0]

# n=3 shell: (level_index, Z_eff², n_quantum, l_quantum)
const _He1_n3_params = [
    (6,  1.236439, 3, 0),
    (7,  1.102898, 3, 0),
    (8,  1.045499, 3, 1),
    (9,  1.001427, 3, 2),
    (10, 0.9926,   3, 1),
]

# Threshold frequency for high-n contributions (n=4)
const _He1_high_n_threshold = 1.25408e16  # Hz

# Inner-shell ionization thresholds
const _He1_ELIM2 = 527490.06   # cm⁻¹
const _He1_ELIM3 = 588451.59   # cm⁻¹

# Group 1: excited state → He II n=2
const _He1_inner_n2 = [(5, 171135.000), (4, 169087.000),
                        (3, 166277.546), (2, 159856.069)]
# Group 2: excited state → He II n=3
const _He1_inner_n3 = [(10, 186209.471), (9, 186101.000), (8, 185564.000),
                        (7, 184864.000),  (6, 183236.000)]


# ============================================================================
# Cross-section tables for the ground state (CROSSHE)
# ============================================================================

const _He1_ground_X505 = Float64[
    7.58, 7.46, 7.33, 7.19, 7.06, 6.94, 6.81,
    6.68, 6.55, 6.43, 6.30, 6.18, 6.05, 5.93,
    5.81, 5.69, 5.57, 5.45, 5.33, 5.21, 5.10,
    4.98, 4.87, 4.76, 4.64, 4.53, 4.42, 4.31,
    4.20, 4.09, 4.00, 3.88, 3.78, 3.68, 3.57,
    3.47, 3.37, 3.27, 3.18, 3.08, 2.98, 2.89,
    2.80, 2.70, 2.61, 2.52, 2.44, 2.35, 2.26,
    2.18, 2.10, 2.02, 1.94, 1.86, 1.78, 1.70,
    1.63, 1.55, 1.48, 1.41, 1.34, 1.28, 1.21,
    1.14, 1.08, 1.02, 0.961, 0.903, 0.847, 0.792,
    0.738, 0.687, 0.637, 0.588, 0.542, 0.497,
    0.454, 0.412, 0.373, 0.335, 0.299, 0.265,
    0.233, 0.202, 0.174, 0.147, 0.123, 0.100,
    0.0795, 0.0609, 0.0443, 0.0315
]

const _He1_ground_X50 = Float64[
    0.0315, 0.0282, 0.0250, 0.0220, 0.0193, 0.0168,
    0.0145, 0.0124, 0.0105, 0.00885, 0.00736,
    0.00604, 0.00489, 0.00389, 0.00303, 0.00231
]

const _He1_ground_X20 = Float64[
    0.00231, 0.00199, 0.00171, 0.00145, 0.00122,
    0.00101, 0.000832, 0.000673, 0.000535, 0.000417,
    0.000318
]

const _He1_ground_X10 = Float64[
    0.000318, 0.000274, 0.000235, 0.000200, 0.000168,
    0.000139, 0.000115, 0.000093, 0.000074, 0.000057,
    0.000044, 0.000032, 0.000023, 0.000016, 0.000010,
    0.000006, 0.000003, 0.000001, 0.0000006,
    0.0000003, 0.0
]

function _He1_ground_bf(ν::Real)
    if ν < 5.945209e15
        return zero(ν)
    end
    λ = _He1_c_ang / ν

    σ_table = if λ > 50.0
        i = clamp(floor(Int, 93.0 - (λ - 50.0) / 5.0), 2, 92)
        λ_node = (92 - i) * 5.0 + 50.0
        frac = (λ - λ_node) / 5.0
        _He1_ground_X505[i] + frac * (_He1_ground_X505[i-1] - _He1_ground_X505[i])
    elseif λ > 20.0
        i = clamp(floor(Int, 17.0 - (λ - 20.0) / 2.0), 2, 16)
        λ_node = (16 - i) * 2.0 + 20.0
        frac = (λ - λ_node) / 2.0
        _He1_ground_X50[i] + frac * (_He1_ground_X50[i-1] - _He1_ground_X50[i])
    elseif λ > 10.0
        i = clamp(floor(Int, 12.0 - (λ - 10.0) / 1.0), 2, 11)
        λ_node = (11 - i) * 1.0 + 10.0
        frac = (λ - λ_node) / 1.0
        _He1_ground_X20[i] + frac * (_He1_ground_X20[i-1] - _He1_ground_X20[i])
    else
        i = clamp(floor(Int, 22.0 - λ / 0.5), 2, 21)
        λ_node = (21 - i) * 0.5
        frac = (λ - λ_node) / 0.5
        _He1_ground_X10[i] + frac * (_He1_ground_X10[i-1] - _He1_ground_X10[i])
    end

    return σ_table * 1.0e-18
end


# ============================================================================
# Cross-section functions for excited states (levels 2–5)
# ============================================================================

# Helper: log-log interpolation for excited-state tables
function _He1_loglog_interp(freq_table::Vector{Float64}, sigma_table::Vector{Float64},
                            ν::Real)
    NP = length(freq_table)
    log_ν = log10(ν)
    i = 2
    while i <= NP && log_ν <= freq_table[i]
        i += 1
    end
    i = clamp(i, 2, NP)
    frac = (log_ν - freq_table[i]) / (freq_table[i-1] - freq_table[i])
    log_σ = sigma_table[i] + frac * (sigma_table[i-1] - sigma_table[i])
    return 10.0^log_σ
end

# --- Level 2: 1s2s ³S ---
const _He1_2s3S_freq = Float64[
    15.956523, 15.923736, 15.888271, 15.849649, 15.807255,
    15.760271, 15.707580, 15.647601, 15.577992, 15.495055,
    15.392451, 15.330345, 15.295609, 15.257851, 15.216496, 15.061770]
const _He1_2s3S_sigma = Float64[
    -18.426022, -18.610700, -18.593051, -18.543304,
    -18.465513, -18.378707, -18.278574, -18.164329,
    -18.033346, -17.882435, -17.705542, -17.605584,
    -17.553459, -17.500667, -17.451318, -17.266686]

function _He1_2s3S_bf(ν::Real)
    ν_threshold = 38454.691 * _He1_c_cgs
    if ν < ν_threshold; return zero(ν); end
    ν_fano = 2.4 * _He1_Rydberg_He * _He1_c_cgs
    if ν > ν_fano
        wavno = ν / _He1_c_cgs
        Ek = (wavno - 38454.691) / _He1_Rydberg_He
        ε = 2.0 * (Ek - 2.47898) / 0.000780
        σ_bg = 0.01521 * (470310.0 / wavno)^3.12 * 8.067e-18
        return σ_bg * (ε - 122.4)^2 / (1.0 + ε^2)
    else
        return _He1_loglog_interp(_He1_2s3S_freq, _He1_2s3S_sigma, ν)
    end
end

# --- Level 3: 1s2s ¹S ---
const _He1_2s1S_freq = Float64[
    15.947182, 15.913654, 15.877320, 15.837666, 15.794025,
    15.745503, 15.690869, 15.628361, 15.555317, 15.467455,
    15.357189, 15.289399, 15.251073, 15.209035, 15.162487, 14.982421]
const _He1_2s1S_sigma = Float64[
    -19.635557, -19.159345, -18.958474, -18.809535,
    -18.676481, -18.546006, -18.410962, -18.264821,
    -18.100205, -17.909165, -17.684370, -17.557867,
    -17.490360, -17.417876, -17.349386, -17.084441]

function _He1_2s1S_bf(ν::Real)
    ν_threshold = 32033.214 * _He1_c_cgs
    if ν < ν_threshold; return zero(ν); end
    ν_fano = 2.4 * _He1_Rydberg_He * _He1_c_cgs
    if ν > ν_fano
        wavno = ν / _He1_c_cgs
        Ek = (wavno - 32033.214) / _He1_Rydberg_He
        ε = 2.0 * (Ek - 2.612316) / 0.00322
        σ_bg = 0.008175 * (484940.0 / wavno)^2.71 * 8.067e-18
        return σ_bg * (ε + 76.21)^2 / (1.0 + ε^2)
    else
        return _He1_loglog_interp(_He1_2s1S_freq, _He1_2s1S_sigma, ν)
    end
end

# --- Level 4: 1s2p ³P ---
const _He1_2p3P_freq = Float64[
    15.943031, 15.909169, 15.872441, 15.832318, 15.788107,
    15.738880, 15.683351, 15.619667, 15.545012, 15.454805,
    15.340813, 15.270195, 15.230054, 15.185821, 15.136567, 14.942557]
const _He1_2p3P_sigma = Float64[
    -19.791021, -19.697886, -19.591421, -19.471855,
    -19.337053, -19.183958, -19.009750, -18.807990,
    -18.570571, -18.288361, -17.943476, -17.738737,
    -17.624154, -17.497163, -17.403183, -17.032999]

function _He1_2p3P_bf(ν::Real)
    ν_threshold = 29223.753 * _He1_c_cgs
    if ν < ν_threshold; return zero(ν); end
    return _He1_loglog_interp(_He1_2p3P_freq, _He1_2p3P_sigma, ν)
end

# --- Level 5: 1s2p ¹P ---
const _He1_2p1P_freq = Float64[
    15.939981, 15.905870, 15.868850, 15.828377, 15.783742,
    15.733988, 15.677787, 15.613218, 15.537343, 15.445346,
    15.328474, 15.255641, 15.214064, 15.168081, 15.116647, 14.911002]
const _He1_2p1P_sigma = Float64[
    -18.798876, -19.685922, -20.011664, -20.143030,
    -20.091354, -19.908333, -19.656788, -19.367745,
    -19.043016, -18.674484, -18.240861, -17.989700,
    -17.852015, -17.702677, -17.525347, -16.816344]

function _He1_2p1P_bf(ν::Real)
    ν_threshold = 27175.76 * _He1_c_cgs
    if ν < ν_threshold; return zero(ν); end
    ν_fano = 2.4 * _He1_Rydberg_He * _He1_c_cgs
    if ν > ν_fano
        wavno = ν / _He1_c_cgs
        Ek = (wavno - 27175.76) / _He1_Rydberg_He
        ε_1S = 2.0 * (Ek - 2.446534) / 0.01037
        ε_1D = 2.0 * (Ek - 2.59427) / 0.00538
        σ_bg = 0.0009487 * (466750.0 / wavno)^3.69 * 8.067e-18
        return σ_bg * ((ε_1S - 29.30)^2 / (1.0 + ε_1S^2) +
                       (ε_1D + 172.4)^2 / (1.0 + ε_1D^2))
    else
        return _He1_loglog_interp(_He1_2p1P_freq, _He1_2p1P_sigma, ν)
    end
end


# ============================================================================
# Hydrogenic cross-section for n=3 shell and high-n levels (Kramers approx)
# ============================================================================

function _He1_hydrogenic_bf(ν::Real, Z_eff2::Real, n::Int, l::Int)
    ν_threshold = Z_eff2 * _He1_Rydberg_He * _He1_c_cgs
    if ν < ν_threshold; return zero(ν); end
    ratio = ν_threshold / ν
    return 7.907e-18 * ratio^3
end


# ============================================================================
# Main function: He I detailed bound-free opacity
# ============================================================================

"""
    He1_detailed_bf!(α, νs, T, number_densities, partition_funcs)

Add detailed He I bound-free opacity to the absorption coefficient vector α.

This is the public interface, designed to be called directly from
`total_continuum_absorption` in ContinuumAbsorption.jl.

Replaces what was previously NO He I bound-free opacity in Korg.
Does NOT include He I free-free, which is already handled by
`positive_ion_ff_absorption!` via `hydrogenic_ff_absorption(ν, T, 2, ...)`.

# Components
  - 5 resolved low-lying states (1s², 1s2s ¹·³S, 1s2p ¹·³P) with dedicated
    cross-section tables and Fano resonance profiles
  - 5 n=3 shell states with fitted effective charges
  - High-n levels (n=4–27) with hydrogenic cross-sections (Z_eff² = 4 - 3/n²)
  - Inner-shell ionization (two-electron transitions)
  - Dissolved-level contribution near the He I series limit

# Arguments
  - `α`: absorption coefficient vector (cm⁻¹), modified in-place
  - `νs`: frequency vector (Hz), sorted
  - `T`: temperature (K)
  - `number_densities`: Dict mapping Species → number density (cm⁻³)
  - `partition_funcs`: Dict mapping Species → partition function callable

# References
  Ported from SYNTHE/ATLAS12 HE1OP (Kurucz 1970, SAO Special Report 309).
  Cross-sections: Mathisen (1984), Fernley et al. (1987), Karzas & Latter (1961).
"""
function He1_detailed_bf!(α::AbstractVector, νs::AbstractVector,
                           T::Real, number_densities::Dict,
                           partition_funcs::Dict)

    nHe_I = get(number_densities, species"He_I", 0.0)
    if nHe_I < 1.0e-30 || T < 100.0
        return α
    end

    # Partition function — needed to convert total n(He I) to level populations
    UHe_I = partition_funcs[species"He_I"](log(T))

    kT_eV = _He1_k_eV * T

    # ──────────────────────────────────────────────────────────────────
    # 1. Boltzmann populations for 10 resolved levels
    #    pop_i = g_i × exp(-χ_i / kT) × n(He I) / U(He I)
    # ──────────────────────────────────────────────────────────────────
    pop = Vector{typeof(T)}(undef, 10)
    for i in 1:10
        pop[i] = _He1_level_g[i] * exp(-_He1_level_χ_eV[i] / kT_eV) * nHe_I / UHe_I
    end

    # ──────────────────────────────────────────────────────────────────
    # 2. Boltzmann populations for high-n levels (n=4–27)
    #    E_n = 24.587 × (1 - 1/n²) eV;  g_n = 4n²
    # ──────────────────────────────────────────────────────────────────
    pop_n = Vector{typeof(T)}(undef, 27)
    for n in 4:27
        E_n = _He1_ionization_eV * (1.0 - 1.0 / n^2)
        g_n = 4.0 * n^2
        pop_n[n] = g_n * exp(-E_n / kT_eV) * nHe_I / UHe_I
    end

    # ──────────────────────────────────────────────────────────────────
    # 3. Dissolved-level populations (BOLTEX and EXLIM)
    #    Smooth opacity from levels between the dissolved threshold
    #    (23.73 eV) and the ionization limit (24.587 eV)
    # ──────────────────────────────────────────────────────────────────
    XR = (nHe_I / UHe_I) * (4.0 / 2.0 / 13.595) * kT_eV
    boltex = exp(-_He1_dissolved_eV / kT_eV) * XR
    exlim  = exp(-_He1_ionization_eV / kT_eV) * XR

    # Precompute the Rydberg frequency for inner-shell threshold calculations
    RYD_freq = _He1_Rydberg_He * _He1_c_cgs

    # ──────────────────────────────────────────────────────────────────
    # 4. Loop over frequencies
    # ──────────────────────────────────────────────────────────────────
    for i in eachindex(νs)
        ν = νs[i]
        freq3 = ν^3

        # Stimulated emission correction
        stim = 1.0 - exp(-_He1_h_eV * ν / kT_eV)

        # --- 4a. Cross-sections for 10 resolved levels ---
        σ1  = ν >= _He1_level_ν_threshold[1]  ? _He1_ground_bf(ν) : zero(T)
        σ2  = ν >= _He1_level_ν_threshold[2]  ? _He1_2s3S_bf(ν)  : zero(T)
        σ3  = ν >= _He1_level_ν_threshold[3]  ? _He1_2s1S_bf(ν)  : zero(T)
        σ4  = ν >= _He1_level_ν_threshold[4]  ? _He1_2p3P_bf(ν)  : zero(T)
        σ5  = ν >= _He1_level_ν_threshold[5]  ? _He1_2p1P_bf(ν)  : zero(T)

        # n=3 shell with fitted Z_eff²
        σ6  = ν >= _He1_level_ν_threshold[6]  ? _He1_hydrogenic_bf(ν, 1.236439, 3, 0) : zero(T)
        σ7  = ν >= _He1_level_ν_threshold[7]  ? _He1_hydrogenic_bf(ν, 1.102898, 3, 0) : zero(T)
        σ8  = ν >= _He1_level_ν_threshold[8]  ? _He1_hydrogenic_bf(ν, 1.045499, 3, 1) : zero(T)
        σ9  = ν >= _He1_level_ν_threshold[9]  ? _He1_hydrogenic_bf(ν, 1.001427, 3, 2) : zero(T)
        σ10 = ν >= _He1_level_ν_threshold[10] ? _He1_hydrogenic_bf(ν, 0.9926,   3, 1) : zero(T)

        # --- 4b. Inner-shell ionization ---
        # Group 1: excited state → He II n=2
        for (lvl_idx, E_level) in _He1_inner_n2
            ν_inner = (_He1_ELIM2 - E_level) * _He1_c_cgs
            if ν >= ν_inner
                Zeff2 = ν_inner / RYD_freq
                σ_inner = _He1_hydrogenic_bf(ν, Zeff2, 1, 0)
                if     lvl_idx == 2; σ2  += σ_inner
                elseif lvl_idx == 3; σ3  += σ_inner
                elseif lvl_idx == 4; σ4  += σ_inner
                elseif lvl_idx == 5; σ5  += σ_inner
                end
            end
        end
        # Group 2: excited state → He II n=3
        for (lvl_idx, E_level) in _He1_inner_n3
            ν_inner = (_He1_ELIM3 - E_level) * _He1_c_cgs
            if ν >= ν_inner
                Zeff2 = ν_inner / RYD_freq
                σ_inner = _He1_hydrogenic_bf(ν, Zeff2, 1, 0)
                if     lvl_idx == 6;  σ6  += σ_inner
                elseif lvl_idx == 7;  σ7  += σ_inner
                elseif lvl_idx == 8;  σ8  += σ_inner
                elseif lvl_idx == 9;  σ9  += σ_inner
                elseif lvl_idx == 10; σ10 += σ_inner
                end
            end
        end

        # --- 4c. Sum bound-free from resolved levels ---
        he1_bf = σ1*pop[1] + σ2*pop[2] + σ3*pop[3] + σ4*pop[4] + σ5*pop[5] +
                 σ6*pop[6] + σ7*pop[7] + σ8*pop[8] + σ9*pop[9] + σ10*pop[10]

        # --- 4d. High-n levels (n=4–27) ---
        if ν >= _He1_high_n_threshold
            for n in 4:27
                Zeff2 = 4.0 - 3.0 / n^2
                σ_n = _He1_hydrogenic_bf(ν, Zeff2, 1, 0)
                he1_bf += σ_n * pop_n[n]
            end
        end

        # --- 4e. Dissolved-level contribution ---
        ex = boltex
        if ν < 2.055e14
            ehvkt = exp(-_He1_h_cgs * ν / (_He1_k_cgs * T))
            ex = exlim / ehvkt
        end
        he1_dissolved = (ex - exlim) * _He1_dissolved_coeff / freq3

        # --- 4f. Total He I bound-free opacity at this frequency ---
        # (no free-free here — that's in positive_ion_ff_absorption!)
        α[i] += (he1_bf + he1_dissolved) * stim
    end

    return α
end