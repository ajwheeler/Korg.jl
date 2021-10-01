module ContinuumAbsorption

using ..Korg: ionization_energies, literals # not sure that this is the best idea
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
    total_continuum_absorption(νs, T, nₑ, number_densities, partition_funcs)

The total continuum linear absoprtion coefficient, α, at many frequencies, ν.

- `νs` are frequencies in Hz
- `T` is temperature in K
- `nₑ` is the electron number density in cm^-3
- `number_densities` is a `Dict` mapping each species to its number density
- `partition_funcs` is a `Dict` mapping each species to its partition function (e.g.
  `Korg.partition_funcs`)
"""
function total_continuum_absorption(νs::AbstractVector{F}, T::F, nₑ::F, number_densities::Dict,
                                    partition_funcs::Dict) where F <: Real
    α = zeros(promote_type(F, valtype(number_densities)), length(νs))

    #TODO check all arguments

    #Hydrogen continuum absorption
    nH_I = number_densities[literals.H_I]
    nH_I_div_U = nH_I / partition_funcs[literals.H_I](T)
    H_I_bf(νs, T, nH_I_div_U; out_α = α)
    H_I_ff(νs, T, number_densities[literals.H_II], nₑ; out_α = α)
    Hminus_bf(νs, T, nH_I_div_U, nₑ; out_α = α)
    Hminus_ff(νs, T, nH_I_div_U, nₑ; out_α = α)
    H2plus_bf_and_ff(νs, T, nH_I_div_U, number_densities[literals.H_II]; out_α = α)

    #He continuum opacities
    He_II_bf(νs, T, number_densities[literals.H_II] / partition_funcs[literals.H_II](T); out_α = α)
    He_II_ff(νs, T, number_densities[literals.He_III], nₑ; out_α = α)
    Heminus_ff(νs, T, number_densities[literals.He_I] / partition_funcs[literals.He_I](T), nₑ;
               out_α = α)

    #electron scattering
    α .+= ContinuumAbsorption.electron_scattering(nₑ)

    α
end


end
