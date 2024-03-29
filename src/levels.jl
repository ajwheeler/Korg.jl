#TODO wrap in module
using HDF5
import Base: ≈, broadcastable
using ProgressMeter

struct Configuration
    config::String
end
broadcastable(e::Configuration) = Ref(e)

# electron configurations are approx equal if one is a terminal substring of the other
function ≈(a::Configuration, b::Configuration) 
    l = min(length(a.config), length(b.config))
    a.config[end-l+1:end] == b.config[end-l+1:end]
end

struct Term
    term::String
end

broadcastable(t::Term) = Ref(t)

function simplify(t::Term)
    t = t.term
    if (length(t) > 0) && islowercase(t[1])
        t = t[2:end]
    end
    if (length(t) > 0) && isdigit(t[end])
        t = t[1:end-1]
    end
    Term(t)
end

function ≈(a::Term, b::Term)
    simplify(a) == simplify(b)
end

struct AtomicLevel
    species::Species
    config::Configuration
    term::Term
    g::UInt16 # 2J + 1
    energy::Float64

    function AtomicLevel(s::Species, config::AbstractString, term::AbstractString, g::Integer, energy::Real)
        new(s, Configuration(config), Term(term), g, energy)
    end
    function AtomicLevel(s::Species, config::Configuration, term::Term, g::Integer, energy::Real)
        new(s, config, term, g, energy)
    end
end

broadcastable(l::AtomicLevel) = Ref(l)

function Base.show(io::IO, level::AtomicLevel)
    print(io, "<")
    show(io, level.species)
    print(io, " ", level.config.config, " ", level.term.term, " g=", level.g, " E=", level.energy, "eV>")
end

"""
TODO
"""
function match_levels(linelist, linelist_levels, coefs; energy_difference_threshold=0.05)
    model_atom_levels, model_atom_bs = coefs

    unique_line_levels = unique([first.(linelist_levels) ; last.(linelist_levels)])
    filter!(unique_line_levels) do level
        !isnothing(level)
    end

    index_vals = @showprogress "matching unique levels" map(unique_line_levels) do level
        match_level(level, model_atom_levels)
    end
    # maps level to (model_atom_ind, level_ind) 
    indices = Dict(l=>i for (l, i) in zip(unique_line_levels, index_vals) if !isnothing(i))

    # maps linelist level to model atom level (for checking how well the matching worked)
    level_matches = map(unique_line_levels) do level
        if level in keys(indices)
            level => model_atom_levels[indices[level]]
        else
            level => nothing
        end
    end |> Dict

    # containers to get filled and returned
    LTE_lines = eltype(linelist)[]
    NLTE_lines = Tuple{Int, Int, eltype(linelist)}[]
    for (line, (lower, upper)) in zip(linelist, linelist_levels)
        if lower in keys(indices) && upper in keys(indices)
            push!(NLTE_lines, (indices[lower], indices[upper], line))
        else
            push!(LTE_lines, line)
        end
    end

    LTE_lines, (NLTE_lines, model_atom_bs), level_matches
end

function match_level(l, model_levels; energy_difference_threshold=0.05)
    species_match = [l.species == m.species for m in model_levels]
    configs_match = [l.config ≈ m.config for m in model_levels]
    terms_match = [l.term == m.term for m in model_levels]
    terms_approx_match = [l.term ≈ m.term for m in model_levels]
    gs_match = [l.g == m.g for m in model_levels]

    for matches in [species_match .& configs_match .& terms_match .& gs_match,
                    species_match .& configs_match .& terms_approx_match .& gs_match,
                    species_match .& configs_match .& terms_match,
                    species_match .& configs_match .& terms_approx_match,
                    species_match .& terms_match,
                    species_match .& terms_approx_match ]
        if sum(matches) == 1
            return findfirst(matches)
        elseif sum(matches) > 1
            ΔE = [abs(m.energy - l.energy) for m in model_levels[matches]]
            # TODO use energy threshold
            #@warn "$level is may correspond to any of $(sum(matches)) levels in the model atom.  Choosing closest_level (ΔE = $(minimum(ΔE)))"
            return findall(matches)[findfirst(ΔE .== minimum(ΔE))]
        end
    end
    nothing
end


"""
The A_X vector for a given [Fe/H] value in 
[the MARCS "standard composition" grid](https://marcs.astro.uu.se/).

Used in [`interpolate_departure_coefs_and_atm`](@ref).
"""
function marcs_standard_composition(Fe_H)
    α_H = clamp(-0.4 * Fe_H, 0, 0.4) + Fe_H
    format_A_X(Fe_H, α_H, Dict("C"=>Fe_H))
end

"""
TODO

TODO default argument vals are an unsatisfactory way to handle this
"""
function interpolate_departure_coefs_and_atm(filenames, Teff=5000.0, logg=4.5, m_H=0.0)
    A_X = marcs_standard_composition(m_H)
    atm = interpolate_marcs(Teff, logg, A_X)

    coefs = map(filenames) do filename
        h5open(filename) do f
            Teffs = read(f, "Teff")
            loggs = read(f, "logg")
            m_Hs = read(f, "m_H")
            x_Fes = read(f, "x_Fe")

            g = Int.(2 * read(f, "J") .+ 1)
            g[g .== -1] .= 0 # handle cases where J was not resolved in model atom 

            levels = AtomicLevel.(Species.(read(f, "species")),
                                  read(f, "configuration"), 
                                  read(f, "term"),
                                  g,
                                  read(f, "E"))

            bs = if occursin("Fe", filename) #TODO
                lazy_multilinear_interpolation(
                    [Teff, logg, m_H],
                    [Teffs, loggs, m_Hs],
                    view(HDF5.readmmap(f["b_array"]), :, :, :, :, :, 1)
                )
            else
                X_Fe = 0 #TODO
                lazy_multilinear_interpolation(
                    [Teff, logg, m_H, X_Fe],
                    [Teffs, loggs, m_Hs, x_Fes],
                    HDF5.readmmap(f["b_array"])
                )
            end

            levels, bs
        end
    end

    coefs = (vcat(first.(coefs)...), hcat(last.(coefs)...))

    atm, coefs
end

"""
TODO
"""
function extract_vald_levels(f)
    lines, _, firstline_ind, shortformat = preprocess_vald_file(f)

    if shortformat
        throw(ArgumentError("Can't extract detailed level info from short-format VALD linelist"))
    end

    @views body = lines[firstline_ind : 4 : end]
    nlines = findfirst(l->l[1]!='\"' || !isuppercase(l[2]), body)-1
    @views body = body[1:nlines]

    # The header is the same for both extract all and extract stellar in "long format"
    CSVheader = ["species", "wl", "loggf", "E_low", "J_lo", "E_up", "J_up", "lower_lande",
                 "upper_lande", "mean_lande", "gamma_rad", "gamma_stark", "gamma_vdW"]
    body = CSV.File(reduce(vcat, codeunits.(body.*"\n")), header=CSVheader, delim=',', 
                    silencewarnings=true)
    
    # the line quantities we need as a vector of named tuples
    transitions = map(eachrow(body)) do row
        (species=Species(row[1].species), E_low=row[1].E_low, J_lo=row[1].J_lo, E_up=row[1].E_up, J_up=row[1].J_up)
    end
    @assert length(transitions) == nlines
                
    level_pairs = map(transitions, 1:nlines) do transition, transition_ind
        line_ind = 4(transition_ind - 1) + firstline_ind 
        lower_level_ind = line_ind + 1
        upper_level_ind = line_ind + 2

        toks_low = split(lines[lower_level_ind][2:end-1])
        @assert length(toks_low) == 3
        lower_level = AtomicLevel(transition.species, toks_low[2], toks_low[3], Int(transition.J_lo*2+1), transition.E_low)

        toks_up = split(lines[upper_level_ind][2:end-1])
        @assert length(toks_up) == 3
        upper_level = AtomicLevel(transition.species, toks_up[2], toks_up[3], Int(transition.J_up*2+1), transition.E_up)

        lower_level, upper_level
    end

    # filter out lines that are filtered out in the linelist
    mask = line_inclusion_criterion.(transitions)

    level_pairs[mask]
end