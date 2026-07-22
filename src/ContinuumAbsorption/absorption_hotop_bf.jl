#=
    absorption_hotop_bf.jl

    Bound-free opacity from doubly and higher ionized metals, ported from
    SYNTHE/ATLAS12 HOTOP.

    This adds 60 photoionization edges for species including C II–V,
    N II–V, O II–V, Ne I–V, Mg II–V, Si II–V, S II–V, and Fe II–V,
    using a modified Seaton formula with parameters from hotop.dat.

    These edges are important for hot stars (T_eff > ~10,000 K) where
    doubly and triply ionized species are abundant. Without this opacity,
    Korg cannot reliably synthesize O/B/early-A star spectra.

    NOTE: The HOTOP *free-free* component (Coulomb Gaunt factors for
    metal ions stages I–V) is already handled by positive_ion_ff_absorption!
    in Korg via Peach (1970) departure coefficients. This module adds
    only the *bound-free* edges.

    == Cross-section formula ==

    The modified Seaton formula for each edge is:
        σ(ν) = σ₀ × (s + ν₀/ν − s·ν₀/ν) × √((ν₀/ν)ⁿ)

    where ν₀ is the threshold frequency, σ₀ is the threshold cross-section,
    s is a shape parameter, and n is a power-law index.

    The opacity contribution from each edge at each atmospheric layer is:
        α_i = σ_i(ν) × g_i × n(species_i) / U(species_i) × exp(−E_exc,i / kT)
              × (1 − exp(−hν/kT))

    == Data file format (hotop.dat) ==

    60 lines, each with 7 parameters:
        ν₀ (Hz), σ₀ (cm²), s, n_power, g_stat, E_exc (eV), species_id

    Species IDs follow the Kurucz numbering convention. The mapping to
    Korg species is provided below.

    == References ==
        Peach, G. 1970, Mem. RAS, 73, 1
        Seaton, M.J. 1958, MNRAS, 118, 504 (original Seaton formula)
        Kurucz, R.L. 1970, SAO Special Report 309
=#

using DelimitedFiles


# ============================================================================
# Kurucz species ID → Korg Species mapping
# ============================================================================

#=
    The Kurucz species numbering for the relevant elements:
    
    Element   Neutral  II    III   IV    V     VI
    ───────   ───────  ──    ───   ──    ─     ──
    C  (Z=6)    21     22    23    24    25    26
    N  (Z=7)    28     29    30    31    32    33
    O  (Z=8)    36     37    38    39    40    41
    Ne (Z=10)   55     56    57    58    59    60
    Mg (Z=12)   78     79    80    81    82    83
    Si (Z=14)  105    106   107   108   109   110
    S  (Z=16)  136    137   138   139   140   141
    Fe (Z=26)  351    352   353   354   355   356
    
    The species_id in hotop.dat refers to the ion whose population
    is used (the species being ionized). For example, species_id = 22
    means C II (singly ionized carbon) is being photoionized → C III.
    
    In Korg's convention, this population is accessed as
    number_densities[species"C_II"] / partition_funcs[species"C_II"](log(T)).
=#

const _hotop_id_to_korg_species = Dict{Int, String}(
    # Carbon
    22 => "C_II",    23 => "C_III",   24 => "C_IV",    25 => "C_V",
    # Nitrogen
    29 => "N_II",    30 => "N_III",   31 => "N_IV",    32 => "N_V",
    # Oxygen
    37 => "O_II",    38 => "O_III",   39 => "O_IV",    40 => "O_V",
    # Neon
    55 => "Ne_I",    56 => "Ne_II",   57 => "Ne_III",  58 => "Ne_IV",
    59 => "Ne_V",    60 => "Ne_VI",
    # Magnesium
    79 => "Mg_II",   80 => "Mg_III",
    # Silicon
    106 => "Si_II",  107 => "Si_III", 108 => "Si_IV",
    # Sulfur
    137 => "S_II",   138 => "S_III",  139 => "S_IV",   140 => "S_V",
    # Iron
    352 => "Fe_II",  353 => "Fe_III", 354 => "Fe_IV",  355 => "Fe_V",
)


# ============================================================================
# Data structures
# ============================================================================

"""
Struct holding the parameters for a single HOTOP bound-free edge.
"""
struct HotopEdge
    ν₀::Float64       # threshold frequency (Hz)
    σ₀::Float64       # threshold cross-section (cm²)
    s::Float64         # shape parameter
    n_power::Int       # power-law exponent (integer)
    g_stat::Float64    # statistical weight of the lower level
    E_exc::Float64     # excitation energy of the lower level (eV)
    species_id::Int    # Kurucz species index
end


# ============================================================================
# Data loading
# ============================================================================

"""
    load_hotop_edges(datadir)

Load the 60 HOTOP bound-free edge parameters from hotop.dat.
Returns a Vector{HotopEdge}.
"""
function load_hotop_edges(datadir::String)
    filepath = joinpath(datadir, "hotop.dat")
    edges = HotopEdge[]

    open(filepath) do f
        for line in eachline(f)
            stripped = lstrip(line)
            if startswith(stripped, '#') || isempty(stripped)
                continue
            end
            vals = parse.(Float64, split(stripped))
            push!(edges, HotopEdge(
                vals[1],           # ν₀
                vals[2],           # σ₀
                vals[3],           # s
                round(Int, vals[4]), # n_power (integer)
                vals[5],           # g_stat
                vals[6],           # E_exc
                round(Int, vals[7]) # species_id
            ))
        end
    end

    return edges
end

const _hotop_edges = load_hotop_edges(_data_dir)

# ============================================================================
# Cross-section evaluation
# ============================================================================

"""
    _hotop_cross_section(edge, ν)

Evaluate the modified Seaton cross-section for a HOTOP edge at frequency ν.
Returns σ in cm². Returns 0 for ν < ν₀ (below threshold).

The formula is:
    σ(ν) = σ₀ × (s + ν₀/ν − s·ν₀/ν) × √((ν₀/ν)^n)
"""
function _hotop_cross_section(edge::HotopEdge, ν::Real)
    if ν < edge.ν₀
        return zero(ν)
    end
    fratio = edge.ν₀ / ν
    σ = edge.σ₀ * (edge.s + fratio - edge.s * fratio) * sqrt(fratio^edge.n_power)
    return σ
end


# ============================================================================
# Main function
# ============================================================================

"""
    hotop_bf_absorption!(α, νs, T, number_densities, partition_funcs)

Add the HOTOP bound-free opacity from doubly+ ionized metals to the
absorption coefficient vector α (cm⁻¹).

This implements the 60 modified-Seaton photoionization edges from
SYNTHE's HOTOP routine, covering C, N, O, Ne, Mg, Si, S, and Fe in
ionization stages II through VI.

# Arguments
  - `α`: absorption coefficient vector (cm⁻¹), modified in-place
  - `νs`: frequency vector (Hz), sorted
  - `T`: temperature (K)
  - `number_densities`: Dict mapping Species → number density (cm⁻³)
  - `partition_funcs`: Dict mapping Species → partition function callable

# Notes
  The SYNTHE HOTOP routine includes a 1% threshold optimization: edges
  contributing less than 1% of the running total are skipped. This
  implementation does NOT include that optimization for simplicity and
  AD compatibility — all edges above threshold are evaluated. The cost
  is negligible (60 edges × simple arithmetic).

  The free-free part of HOTOP (Coulomb Gaunt factors for metal ions
  stages I–V) is already handled by positive_ion_ff_absorption! in Korg.

# References
  - Peach, G. 1970, Mem. RAS, 73, 1
  - Seaton, M.J. 1958, MNRAS, 118, 504
  - Kurucz, R.L. 1970, SAO Special Report 309
"""
function hotop_bf_absorption!(α::AbstractVector, νs::AbstractVector,
                               T::Real, number_densities::Dict,
                               partition_funcs::Dict)

    edges = _hotop_edges

    kT_eV = 8.617333262e-5 * T
    h_eV  = 4.135667696e-15

    # Precompute population factors for each unique species in the edge list
    # pop_factor = n(species) / U(species)
    # Cache to avoid recomputing for edges that share the same species
    pop_cache = Dict{Int, typeof(T)}()

    for edge in edges
        if haskey(pop_cache, edge.species_id)
            continue
        end

        spec_str = get(_hotop_id_to_korg_species, edge.species_id, nothing)
        if spec_str === nothing
            pop_cache[edge.species_id] = zero(T)
            continue
        end

        spec = Species(spec_str)
        ndens = get(number_densities, spec, 0.0)
        if ndens <= 0.0
            pop_cache[edge.species_id] = zero(T)
            continue
        end

        U = partition_funcs[spec](log(T))
        pop_cache[edge.species_id] = ndens / U
    end

    # Loop over frequencies
    for i in eachindex(νs)
        ν = νs[i]

        # Stimulated emission correction
        stim = 1.0 - exp(-h_eV * ν / kT_eV)

        # Sum contributions from all edges above threshold
        α_hotop = zero(T)
        for edge in edges
            if ν < edge.ν₀
                continue
            end

            pop = pop_cache[edge.species_id]
            if pop <= 0.0
                continue
            end

            σ = _hotop_cross_section(edge, ν)
            # Boltzmann factor for excitation + statistical weight
            boltz = edge.g_stat * exp(-edge.E_exc / kT_eV)
            α_hotop += σ * pop * boltz
        end

        α[i] += α_hotop * stim
    end

    return α
end