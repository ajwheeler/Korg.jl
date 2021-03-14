module ContinuumOpacity

using ..SSSynth: ionization_energies # not sure that this is the best idea
include("../constants.jl") # I'm not thrilled to duplicate this, but I think it's probably alright

# define helper functions
include("hydrogenic_bf_ff.jl")

export H_I_bf, H_I_ff, Hminus_bf, Hminus_ff, H2plus_bf_and_ff
include("opacity_H.jl")

export He_II_bf, He_II_ff, Heminus_ff
include("opacity_He.jl")

end
