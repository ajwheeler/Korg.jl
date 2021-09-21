module ContinuumOpacity

using ..Korg: ionization_energies # not sure that this is the best idea
include("../constants.jl") # I'm not thrilled to duplicate this, but I think it's probably alright

# define helper functions
include("hydrogenic_bf_ff.jl")

export H_I_bf, H_I_ff, Hminus_bf, Hminus_ff, H2plus_bf_and_ff
include("opacity_H.jl")

export He_II_bf, He_II_ff, Heminus_ff
include("opacity_He.jl")


# the following are only imported for computing experimental metal bf continuum opacities
using ..Korg: partition_funcs, Species, Formula, ismolecule, get_roman_numerals
export absorption_coef_bf_TOPBase
include("opacity_metal.jl")


"""
    electron_scattering(nₑ, ρ)

Compute the opacity from scattering off of free electrons. This has no wavelength dependence. It
assumes isotropic scattering

# Arguments
- `nₑ::F`: number density of free electrons (in cgs)
- `ρ::F`: the mass density

# Notes
I adopted the formula described in section 5.12 of Kurucz (1970) and the equation in the electron
scattering subsection of Gray (2005); the actual coefficient value comes from the latter. It turns
out that the coefficient in Kurucz (1970) has a typo (it's a factor of 10 too large).
"""
electron_scattering(nₑ::F, ρ::F) where {F<:Real} = 0.6648e-24*nₑ/ρ
end
