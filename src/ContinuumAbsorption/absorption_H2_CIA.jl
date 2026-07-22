#=
    absorption_H2_CIA.jl

    H₂ collision-induced absorption (CIA) opacity from H₂-H₂ and H₂-He
    collisions.

    Ported from SYNTHE/ATLAS12 H2COLLOP subroutine using tabulated absorption
    coefficients from Borysow, Jorgensen & Zheng (1997, A&A 324, 185).

    Two 2D tables are read from h2collop.dat:
        H2HE(7, 81): log₁₀ of H₂-He CIA absorption coefficient
        H2H2(7, 81): log₁₀ of H₂-H₂ CIA absorption coefficient

    Grid axes:
        Temperature: T = 1000, 2000, ..., 7000 K  (7 points, 1000 K steps)
        Wavenumber:  ν̃ = 0, 250, ..., 20000 cm⁻¹  (81 points, 250 cm⁻¹ steps)

    The final opacity is:
        α_CIA = (10^σ_H2He · n(He I) + 10^σ_H2H2 · n(H₂)) · n(H₂) · (1 - exp(-hν/kT))

    where σ values are bilinearly interpolated from the tables.

    NOTE: The original Kurucz Fortran code had a bug in the temperature
    interpolation weights (swapped low/high), which was corrected in the F90
    modernization. This Julia implementation uses the corrected weights.

    References:
        Borysow, A., Jorgensen, U.G. & Zheng, C. 1997, A&A 324, 185
        Borysow, A. 2002, A&A 390, 779 (updated H₂-H₂ tables)
=#

using DelimitedFiles

# ============================================================================
# Data loading
# ============================================================================

"""
    _load_h2_CIA_tables(datadir)

Load the H₂-He and H₂-H₂ CIA coefficient tables from h2collop.dat.
Returns a tuple (H2HE, H2H2) where each is a 7×81 matrix of log₁₀(σ_CIA).

The file format is:
  - Header lines starting with '#'
  - 81 rows of 7 values for H₂-He (temperatures across columns, wavenumbers down rows)
  - 81 rows of 7 values for H₂-H₂
"""
function _load_h2_CIA_tables(datadir::String)
    filepath = joinpath(datadir, "h2collop.dat")

    # Read all non-comment lines
    lines = String[]
    open(filepath) do f
        for line in eachline(f)
            stripped = lstrip(line)
            if !startswith(stripped, '#') && !isempty(stripped)
                push!(lines, stripped)
            end
        end
    end

    # Parse into a single matrix: should be 162 rows × 7 columns
    # First 81 rows = H2HE, next 81 rows = H2H2
    all_data = zeros(162, 7)
    for (i, line) in enumerate(lines)
        vals = parse.(Float64, split(line))
        all_data[i, :] = vals
    end

    # Transpose to get (temperature_index, wavenumber_index) = (7, 81)
    H2HE = Matrix{Float64}(all_data[1:81, :]')     # 7×81
    H2H2 = Matrix{Float64}(all_data[82:162, :]')    # 7×81

    return H2HE, H2H2
end

const H2HE_table, H2H2_table = _load_h2_CIA_tables(_data_dir)

# ============================================================================
# Interpolation
# ============================================================================

"""
    _interp_H2_CIA(table, T, wavenumber)

Bilinearly interpolate a CIA coefficient table at temperature T (K)
and wavenumber ν̃ (cm⁻¹). Returns the CIA coefficient (NOT log₁₀).

The interpolation is done in log₁₀ space (linear interpolation of log₁₀ values),
then exponentiated.

Temperature grid: 1000, 2000, ..., 7000 K (7 points)
Wavenumber grid: 0, 250, ..., 20000 cm⁻¹ (81 points)

Returns 0.0 for wavenumbers > 20,000 cm⁻¹ or temperatures outside 1000–7000 K.
"""
function _interp_H2_CIA(table::Matrix{Float64}, T::Real, wavenumber::Real)
    # Out of range checks
    if wavenumber > 20000.0 || wavenumber < 0.0 || T < 1000.0 || T > 7000.0
        return zero(promote_type(typeof(T), typeof(wavenumber)))
    end

    # Wavenumber index: ν̃ = 0, 250, ..., 20000 → index 1..81
    nu_idx = wavenumber / 250.0 + 1.0
    nu_lo = clamp(floor(Int, nu_idx), 1, 80)
    frac_nu = nu_idx - nu_lo

    # Temperature index: T = 1000, 2000, ..., 7000 → index 1..7
    t_idx = (T - 1000.0) / 1000.0 + 1.0
    t_lo = clamp(floor(Int, t_idx), 1, 6)
    frac_t = t_idx - t_lo

    # Bilinear interpolation in log₁₀ space
    # First interpolate in wavenumber at both temperature bracketing points
    log_val_t_lo = table[t_lo, nu_lo] * (1.0 - frac_nu) + table[t_lo, nu_lo+1] * frac_nu
    log_val_t_hi = table[t_lo+1, nu_lo] * (1.0 - frac_nu) + table[t_lo+1, nu_lo+1] * frac_nu

    # Then interpolate in temperature
    log_val = log_val_t_lo * (1.0 - frac_t) + log_val_t_hi * frac_t

    return 10.0^log_val
end


# ============================================================================
# Public interface
# ============================================================================

"""
    H2_CIA_absorption!(α, νs, T, number_densities)

Add the H₂ collision-induced absorption (CIA) opacity to the absorption
coefficient vector α (cm⁻¹).

CIA arises from transient electric dipole moments induced during H₂-H₂ and
H₂-He collisions, producing broadband absorption primarily in the infrared.
This opacity source is important for cool dwarfs (T_eff < 4000 K) and
substellar objects where molecular hydrogen is abundant.

# Arguments
  - `α`: absorption coefficient vector to be modified in-place (cm⁻¹)
  - `νs`: frequency vector (Hz), sorted
  - `T`: temperature (K)
  - `number_densities`: Dict mapping Species → number density (cm⁻³)

# Physics

The CIA absorption coefficient is:
```
α_CIA = (σ_H2He · n(He I) + σ_H2H2 · n(H₂)) · n(H₂) · (1 - exp(-hν/kT))
```
where σ values are interpolated from the Borysow, Jorgensen & Zheng (1997)
tables. Note the n(H₂)² dependence for the H₂-H₂ component — CIA scales
as the product of both collider densities.

# Notes

In SYNTHE's H2COLLOP, the assembly (within COOLOP) is:
```fortran
AH2COLL(J) = (10^XH2HE * XNF(J,3) + 10^XH2H2 * XNH2(J)) * XNH2(J) / RHO(J) * STIM(J)
```
where XNF(J,3) = n(He I)/ρ and XNH2(J) = n(H₂). Here we use Korg's number
densities directly (which give n in cm⁻³) and do not divide by ρ because
Korg's α is in cm⁻¹ (linear absorption coefficient), not cm²/g (mass
absorption coefficient).

The species keys used are `species"H2 I"` for molecular hydrogen and
`species"He I"` for neutral helium. If H₂ is not present in the number
densities dict, no CIA contribution is added.

# References
  - Borysow, A., Jorgensen, U.G. & Zheng, C. 1997, A&A 324, 185
"""
function H2_CIA_absorption!(α::AbstractVector, νs::AbstractVector,
                             T::Real, number_densities::Dict)

    # Physical constants
    h_eV = 4.135667696e-15   # Planck constant in eV·s
    k_eV = 8.617333262e-5    # Boltzmann constant in eV/K
    c_cm = 2.99792458e10     # speed of light in cm/s

    # Get number densities
    nH2 = get(number_densities, species"H2", 0.0)
    if nH2 <= 0.0
        return α
    end
    nHe = get(number_densities, species"He I", 0.0)

    # Temperature range check: tables only cover 1000–7000 K
    if T < 1000.0 || T > 7000.0
        return α
    end

    for i in eachindex(νs)
        # Convert frequency to wavenumber (cm⁻¹)
        wavenumber = νs[i] / c_cm

        # Skip if above table range
        if wavenumber > 20000.0
            continue
        end

        # Interpolate CIA coefficients from tables
        σ_H2He = _interp_H2_CIA(H2HE_table, T, wavenumber)
        σ_H2H2 = _interp_H2_CIA(H2H2_table, T, wavenumber)

        # Stimulated emission correction
        E_eV = h_eV * νs[i]
        stim = 1.0 - exp(-E_eV / (k_eV * T))

        # CIA absorption: note the n(H₂)² dependence for H₂-H₂
        α[i] += (σ_H2He * nHe + σ_H2H2 * nH2) * nH2 * stim
    end

    return α
end