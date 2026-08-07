module ContinuumAbsorption
export total_continuum_absorption

using ..Korg: ionization_energies, Species, @species_str, _data_dir # not sure that this is the best idea
using ..Korg: Interval, closed_interval, contained, contained_slice, λ_to_ν_bound, hummer_mihalas_w
include("../constants.jl") # I'm not thrilled to duplicate this, but I think it's probably alright

include("bounds_checking.jl") # define helper functions
include("hydrogenic_bf_ff.jl")
include("absorption_H.jl")
include("absorption_He.jl")
include("absorption_ff_positive_ion.jl")
include("absorption_metals_bf.jl")
include("scattering.jl")
include("absorption_mol_photodissociation.jl")
include("absorption_H2_CIA.jl")
include("absorption_He1_detailed.jl")
include("absorption_hotop_bf.jl")

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
function total_continuum_absorption(νs, T, nₑ, number_densities::Dict, partition_funcs::Dict;
                                    error_oobounds=false)
    α = zeros(promote_type(eltype(νs), typeof(T), typeof(nₑ), valtype(number_densities)),
              length(νs))

    kwargs = Dict(:out_α => α, :error_oobounds => error_oobounds)

    # used more than once
    nH_I = number_densities[species"H_I"]
    invU_H_I = 1 / partition_funcs[species"H I"](log(T))

    # Hydrogen continuum absorption
    # note: inclusion of He I ndens below is NOT a typo
    α .+= H_I_bf(νs, T, nH_I, number_densities[species"He I"], nₑ, invU_H_I)

    Hminus_bf(νs, T, number_densities[species"H-"], nₑ; kwargs...)
    Hminus_ff(νs, T, nH_I * invU_H_I, nₑ; kwargs...)
    H2plus_bf_and_ff(νs, T, nH_I, number_densities[species"H_II"]; kwargs...)

    # He continuum absorption isn't actually important, but here we are
    Heminus_ff(νs, T, number_densities[species"He_I"] / partition_funcs[species"He_I"](log(T)), nₑ;
               kwargs...)

    # ff absorption where participating species are positive ions 
    # i.e. H I ff is included but not H⁻ ff or He⁻ ff 
    positive_ion_ff_absorption!(α, νs, T, number_densities, nₑ)

    # bf absorption by metals from TOPBase and NORAD
    metal_bf_absorption!(α, νs, T, number_densities)

    # He I detailed bound-free (HE1OP port from SYNTHE)
    # This adds opacity from 10 resolved He I states + high-n levels +
    # inner-shell ionization + dissolved levels near the series limit.
    # He I free-free is NOT included here (it's in positive_ion_ff_absorption!).
    He1_detailed_bf!(α, νs, T, number_densities, partition_funcs)

    # molecular photodissociation (specifically OH and CH)
    mol_photodissociation_absorption!(α, νs, T, number_densities)

    # H2 collision-induced absorption (CIA)
    H2_CIA_absorption!(α, νs, T, number_densities)

    # HOTOP bound-free absorption from doubly+ ionized metals
    hotop_bf_absorption!(α, νs, T, number_densities, partition_funcs)

    # scattering
    α .+= electron_scattering(nₑ)
    α .+= rayleigh(νs, nH_I, number_densities[species"He_I"], number_densities[species"H2"])

    α
end

end
