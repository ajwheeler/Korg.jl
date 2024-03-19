#TODO wrap in module

import Base: ≈, broadcastable

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
    config::Configuration
    term::Term
    g::UInt8 # 2J + 1
    energy::Float64

    function AtomicLevel(config::String, term::String, g::Integer, energy::Real)
        new(Configuration(config), Term(term), g, energy)
    end
    function AtomicLevel(config::Configuration, term::Term, g::Integer, energy::Real)
        new(config, term, g, energy)
    end
end

broadcastable(l::AtomicLevel) = Ref(l)

function Base.show(io::IO, ::MIME"text/plain", level::AtomicLevel)
    print(io, level.config.config, " ", level.term.term, " g=", level.g, " E=", level.energy, "eV")
end

"""
TODO
"""
function match_levels(linelist_levels, model_levels; energy_difference_threshold=0.05)
    map(linelist_levels) do level
        configs_match = [level.config ≈ m.config for m in model_levels]
        terms_match = [level.term == m.term for m in model_levels]
        terms_approx_match = [level.term ≈ m.term for m in model_levels]
        gs_match = [level.g == m.g for m in model_levels]

        for (i, matches) in enumerate([
                configs_match .& terms_match .& gs_match,
                configs_match .& terms_approx_match .& gs_match,
                configs_match .& terms_match,
                configs_match .& terms_approx_match,
                terms_match,
                terms_approx_match
                ])
            
            if sum(matches) == 1
                return i, findfirst(matches)
            elseif sum(matches) > 1
                ΔE = [abs(m.energy - level.energy) for m in model_levels[matches]]
                # TODO use energy threshold
                #@warn "$level is may correspond to any of $(sum(matches)) levels in the model atom.  Choosing closest_level (ΔE = $(minimum(ΔE)))"
                return (i + 0.5, findfirst(ΔE .== minimum(ΔE)))
            end
        end
        0, nothing
    end
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

using HDF5
"""
TODO
"""
function interpolate_departure_coefs_and_atm(Teff, logg, m_H, filenames)
    A_X = marcs_standard_composition(m_H)
    atm = interpolate_marcs(Teff, logg, A_X)
    n_layers = length(atm.layers)

    coefs = map(filenames) do filename
        h5open(filename) do f
            Teffs = read(f, "Teff")
            loggs = read(f, "logg")
            m_Hs = read(f, "m_H")
            x_Fes = read(f, "x_Fe")

            g = Int.(2 * read(f, "J") .+ 1)
            g[g .== -1] .= 0 # handle cases where J was not resolved in model atom 

            levels = AtomicLevel.(read(f, "configuration"), 
                                  read(f, "term"),
                                  g,
                                  read(f, "E"))
            n_levels = length(levels)
            specs = fill(species"Fe I", n_levels) #TODO

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

            specs, levels, bs
        end
    end

    atm, coefs
end