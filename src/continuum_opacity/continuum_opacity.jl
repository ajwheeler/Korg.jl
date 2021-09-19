module ContinuumOpacity

using ..Korg: ionization_energies # not sure that this is the best idea
include("../constants.jl") # I'm not thrilled to duplicate this, but I think it's probably alright

# define helper functions
include("hydrogenic_bf_ff.jl")

export abs_H_I_bf, abs_H_I_ff, abs_Hminus_bf, abs_Hminus_ff, abs_H2plus_bf_and_ff
include("opacity_H.jl")

export abs_He_II_bf, abs_He_II_ff, abs_Heminus_ff
include("opacity_He.jl")

include("scattering.jl")

# the following are only imported for computing experimental metal bf continuum opacities
using ..Korg: partition_funcs, Species, Formula, ismolecule, get_roman_numeral
export absorption_coef_bf_TOPBase
include("opacity_metal.jl")

end
