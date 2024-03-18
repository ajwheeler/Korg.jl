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