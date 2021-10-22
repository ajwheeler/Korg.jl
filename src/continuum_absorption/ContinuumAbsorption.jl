module ContinuumAbsorption

using ..Korg: ionization_energies, literals, _data_dir # not sure that this is the best idea
include("../constants.jl") # I'm not thrilled to duplicate this, but I think it's probably alright

# define helper functions
include("bounds_checking.jl")
include("hydrogenic_bf_ff.jl")

export H_I_bf, H_I_ff, Hminus_bf, Hminus_ff, H2plus_bf_and_ff
include("absorption_H.jl")

export He_II_bf, He_II_ff, Heminus_ff
include("absorption_He.jl")

include("scattering.jl")

# the following are only imported for computing experimental metal bf continuum opacities
using ..Korg: partition_funcs, Species, Formula, ismolecule, get_roman_numeral
export absorption_coef_bf_TOPBase
include("absorption_metal.jl")

export total_continuum_absorption

"""
    total_continuum_absorption(νs, T, nₑ, number_densities, partition_funcs; extapolate_bc)

The total continuum linear absoprtion coefficient, α, at many frequencies, ν.

# Arguments

- `νs` are frequencies in Hz
- `T` is temperature in K
- `nₑ` is the electron number density in cm^-3
- `number_densities` is a `Dict` mapping each species to its number density
- `partition_funcs` is a `Dict` mapping each species to its partition function (e.g.
  `Korg.partition_funcs`)
- `extrapolate_bc` specifies how most continnum absorption sources passed frequencies and 
   temperature values that are out of bounds for their implementation. The default value is `0.0`, 
   which simply ignores an opacity source at those values. For additional details see the bullet 
   point at [Continuum Absorption Kwargs](@ref)

!!! note
    For efficiency reasons, `νs` must be sorted. While this function technically supports any 
    sorted `AbstractVector`, it is most effient when passed an  `AbstractRange`.
"""
function total_continuum_absorption(νs::AbstractVector{F}, T::F, nₑ::F, number_densities::Dict,
                                    partition_funcs::Dict; extrapolate_bc = 0.0) where F <: Real
    α = zeros(promote_type(F, valtype(number_densities)), length(νs))

    #TODO check all arguments

    kwargs = Dict(:out_α=>α, :extrapolate_bc=>extrapolate_bc)

    #Hydrogen continuum absorption
    nH_I = number_densities[literals.H_I]
    nH_I_div_U = nH_I / partition_funcs[literals.H_I](T)
    H_I_bf(νs, T, nH_I_div_U; kwargs...)
    H_I_ff(νs, T, number_densities[literals.H_II], nₑ; kwargs...)
    Hminus_bf(νs, T, nH_I_div_U, nₑ; kwargs...)
    Hminus_ff(νs, T, nH_I_div_U, nₑ; kwargs...)
    H2plus_bf_and_ff(νs, T, nH_I_div_U, number_densities[literals.H_II]; kwargs...)

    #He continuum opacities
    He_II_bf(νs, T, number_densities[literals.H_II] / partition_funcs[literals.H_II](T); kwargs...)
    He_II_ff(νs, T, number_densities[literals.He_III], nₑ; kwargs...)
    Heminus_ff(νs, T, number_densities[literals.He_I] / partition_funcs[literals.He_I](T), nₑ;
               kwargs...)

    #electron scattering
    α .+= ContinuumAbsorption.electron_scattering(nₑ)

    α
end


end
