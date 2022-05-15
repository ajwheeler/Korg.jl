module ContinuumAbsorption

using ..Korg: ionization_energies, @species_str, _data_dir # not sure that this is the best idea
using ..Korg: Interval, closed_interval, contained, contained_slice, λ_to_ν_bound
include("../constants.jl") # I'm not thrilled to duplicate this, but I think it's probably alright

# define helper functions
include("bounds_checking.jl")
include("hydrogenic_bf_ff.jl")

include("absorption_H.jl")
include("absorption_He.jl")

include("absorption_ff_positive_ion.jl")
include("absorption_ff_neutral_metals_molecules.jl")
include("scattering.jl")

# the following are only imported for computing experimental metal bf continuum opacities
using ..Korg: partition_funcs, Species, Formula, ismolecule, get_roman_numeral
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
    nH_I_div_U = nH_I / partition_funcs[species"H_I"](log(T))

    # Hydrogen continuum absorption
    H_I_bf(νs, T, nH_I_div_U; kwargs...)
    Hminus_bf(νs, T, nH_I_div_U, nₑ; kwargs...)
    Hminus_ff(νs, T, nH_I_div_U, nₑ; kwargs...)
    H2plus_bf_and_ff(νs, T, nH_I, number_densities[species"H_II"]; kwargs...)

    # He continuum absorption
    He_II_bf(νs, T, number_densities[species"H_II"]/partition_funcs[species"H_II"](log(T)); kwargs...)
    Heminus_ff(νs, T, number_densities[species"He_I"] / partition_funcs[species"He_I"](log(T)), nₑ;
               kwargs...)

    # ff absorption where participating species are positive ions 
    # i.e. H I ff is included but not H⁻ ff or He⁻ ff 
    positive_ion_ff_absorption!(α, νs, T, number_densities, nₑ)

    # ff absorption where participating species are neutral metals
    Ominus_ff(α, νs, T, get(number_densities, species"O_I", 0.0), nₑ)

    # there is some minor cause for concern with these 2 sources: set the unit test for Nminus_ff
    #Cminus_ff(α, νs, T, get(number_densities, species"C_I", 0.0), nₑ)
    #Nminus_ff(α, νs, T, get(number_densities, species"N_I", 0.0), nₑ)

    # ff absorption where participating species are neutral molecules
    H2minus_ff(α, νs, T, get(number_densities, species"H2_I", 0.0), nₑ)
    COminus_ff(α, νs, T, get(number_densities, species"CO_I", 0.0), nₑ)

    # scattering
    α .+= electron_scattering(nₑ)
    α .+= rayleigh(νs, nH_I, number_densities[species"He_I"], number_densities[species"H2"])

    α
end


end
