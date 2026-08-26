module ContinuumAbsorption
export total_continuum_absorption

using ..Korg: ionization_energies, Species, @species_str, _data_dir # not sure that this is the best idea
using ..Korg: Interval, closed_interval, contained, contained_slice, λ_to_ν_bound, hummer_mihalas_w
include("../constants.jl") # I'm not thrilled to duplicate this, but I think it's probably alright

include("hydrogenic_bf_ff.jl")
include("absorption_H.jl")
include("absorption_He.jl")
include("absorption_ff_positive_ion.jl")
include("absorption_metals_bf.jl")
include("scattering.jl")

"""
    total_continuum_absorption(νs, T, nₑ, number_densities, partition_funcs)

The total continuum linear absoprtion coefficient, α, at many frequencies, ν.

# Arguments

  - `νs` are frequencies in Hz. `νs` must be sorted. While this function
    technically supports any sorted `AbstractVector`, it is most effient when
    passed an  `AbstractRange`.
  - `T` is temperature in K
  - `nₑ` is the electron number density in cm^-3
  - `number_densities` is a `Dict` mapping each `Species` to its number density
  - `partition_funcs` is a `Dict` mapping each `Species` to its partition function (e.g.
    `Korg.partition_funcs`)
"""
function total_continuum_absorption(νs, T, nₑ, number_densities::Dict, partition_funcs::Dict)
    α = zeros(promote_type(eltype(νs), typeof(T), typeof(nₑ), valtype(number_densities)),
              length(νs))

    # used more than once
    nH_I = number_densities[species"H_I"]
    invU_H_I = 1 / partition_funcs[species"H I"](log(T))

    # Hydrogen continuum absorption
    # note: inclusion of He I ndens below is NOT a typo
    α .+= H_I_bf(νs, T, nH_I, number_densities[species"He I"], nₑ, invU_H_I)

    Hminus_bf!(α, νs, T, number_densities[species"H-"], nₑ)
    Hminus_ff!(α, νs, T, nH_I * invU_H_I, nₑ)
    H2plus_bf_and_ff!(α, νs, T, nH_I, number_densities[species"H_II"])

    # He continuum absorption isn't actually important, but here we are
    Heminus_ff!(α, νs, T,
                number_densities[species"He_I"] / partition_funcs[species"He_I"](log(T)), nₑ)

    # ff absorption where participating species are positive ions 
    # i.e. H I ff is included but not H⁻ ff or He⁻ ff 
    positive_ion_ff_absorption!(α, νs, T, number_densities, nₑ)

    # bf absorption by metals from TOPBase and NORAD
    metal_bf_absorption!(α, νs, T, number_densities)

    # scattering
    α .+= electron_scattering(nₑ)
    α .+= rayleigh(νs, nH_I, number_densities[species"He_I"], number_densities[species"H2"])

    α
end

end
