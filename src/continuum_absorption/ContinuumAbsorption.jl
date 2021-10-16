module ContinuumAbsorption

using ..Korg: ionization_energies, @species_str, _data_dir # not sure that this is the best idea
using ..Korg: Interval, closed_interval, contained, contained_slice, λ_to_ν_bound
include("../constants.jl") # I'm not thrilled to duplicate this, but I think it's probably alright

# define helper functions
include("bounds_checking.jl")
include("hydrogenic_bf_ff.jl")
include("stancil_tables.jl")

export H_I_bf, H_I_ff, Hminus_bf, Hminus_ff, H2plus_bf_and_ff
include("absorption_H.jl")

export He_II_bf, He_II_ff, Heminus_ff
include("absorption_He.jl")

include("absorption_ff_metal.jl")
include("scattering.jl")

# the following are only imported for computing experimental metal bf continuum opacities
using ..Korg: partition_funcs, Species, Formula, ismolecule, get_roman_numeral
export absorption_coef_bf_TOPBase
include("absorption_metal.jl")

export total_continuum_absorption

"""
    total_continuum_absorption(νs, T, nₑ, number_densities, partition_funcs; error_oobounds)

The total continuum linear absoprtion coefficient, α, at many frequencies, ν.

# Arguments

- `νs` are frequencies in Hz
- `T` is temperature in K
- `nₑ` is the electron number density in cm^-3
- `number_densities` is a `Dict` mapping each `Species` to its number density
- `partition_funcs` is a `Dict` mapping each `Species` to its partition function (e.g.
  `Korg.partition_funcs`)
- `error_oobounds::Bool` specifies the behavior of most continnum absorption sources when passed
   frequencies or temperature values that are out of bounds for their implementation. When `false`
   (the default), those absorption sources are ignored at those values. Otherwise, an error is
   thrown.

!!! note
    For efficiency reasons, `νs` must be sorted. While this function technically supports any 
    sorted `AbstractVector`, it is most effient when passed an  `AbstractRange`.
"""
function total_continuum_absorption(νs::AbstractVector{F}, T::F, nₑ::F, number_densities::Dict,
                                    partition_funcs::Dict; error_oobounds = false) where F <: Real
    α = zeros(promote_type(F, valtype(number_densities)), length(νs))

    kwargs = Dict(:out_α => α, :error_oobounds => error_oobounds)

    #parameters used more than once
    nH_I = number_densities[species"H_I"]
    nH_I_div_U = nH_I / partition_funcs[species"H_I"](T)

    # Hydrogen continuum absorption
    H_I_bf(νs, T, nH_I_div_U; kwargs...)
    H_I_ff(νs, T, number_densities[species"H_II"], nₑ; kwargs...)
    Hminus_bf(νs, T, nH_I_div_U, nₑ; kwargs...)
    Hminus_ff(νs, T, nH_I_div_U, nₑ; kwargs...)
    H2plus_bf_and_ff(νs, T, nH_I, number_densities[species"H_II"]; kwargs...)

    # He continuum absorption
    He_II_bf(νs, T, number_densities[species"H_II"] / partition_funcs[species"H_II"](T); kwargs...)
    He_II_ff(νs, T, number_densities[species"He_III"], nₑ; kwargs...)
    Heminus_ff(νs, T, number_densities[species"He_I"] / partition_funcs[species"He_I"](T), nₑ;
               kwargs...)

    # scattering
    α .+= electron_scattering(nₑ)
    α .+= rayleigh(νs, nH_I, number_densities[species"He_I"], number_densities[species"H2"])

    α
end


end
