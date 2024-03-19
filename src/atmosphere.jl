# used when downloading model atmosphere archive
using ProgressMeter: Progress, update!, finish!
using Pkg.Artifacts: @artifact_str
 
abstract type ModelAtmosphere end

struct PlanarAtmosphereLayer{F1, F2, F3, F4, F5}
    tau_5000::F1                #dimensionless (used for legacy radiative transfer)
    z::F2                       #cm
    temp::F3                    #K
    electron_number_density::F4 #cm^-3
    number_density::F5          #cm^-3
end

struct PlanarAtmosphere{F1, F2, F3, F4, F5} <: ModelAtmosphere
    layers::Vector{PlanarAtmosphereLayer{F1, F2, F3, F4, F5}}
end

struct ShellAtmosphereLayer{F1, F2, F3, F4, F5}
    tau_5000::F1                #dimensionless (used for legacy radiative transfer)
    z::F2                       #cm
    temp::F3                    #K
    electron_number_density::F4 #cm^-3
    number_density::F5          #cm^-3
end

struct ShellAtmosphere{F1, F2, F3, F4, F5, F6} <: ModelAtmosphere
    layers::Vector{ShellAtmosphereLayer{F1, F2, F3, F4, F5}}
    R::F6 #the radius of the star where τ_ros == 1, i.e. the photosphere (not the top)
end

"""
    PlanarAtmosphere(atm::ShellAtmosphere)

Construct a planar atmosphere with the data from a shell atmosphere.  Mostly useful for testing.
"""
function PlanarAtmosphere(atm::ShellAtmosphere)
    PlanarAtmosphere([PlanarAtmosphereLayer(l.tau_5000, l.z, l.temp, l.electron_number_density, 
                                            l.number_density) for l in atm.layers])
end

"""
    ShellAtmosphere(atm::PlanarAtmosphere, R)

Construct a shell atmosphere with the data from a planar atmosphere and an outer radius.  Mostly 
useful for testing.
"""
function ShellAtmosphere(atm::PlanarAtmosphere, R)
    ShellAtmosphere([ShellAtmosphereLayer(l.tau_5000, l.z, l.temp, l.electron_number_density, 
                                          l.number_density) for l in atm.layers], R)
end

#pretty-printing
function Base.show(io::IO, m::MIME"text/plain", atm::A) where A <: ModelAtmosphere
    print(io, "$(A) with $(length(atm.layers)) layers")
end

"""   
    tau_5000s(atm::ModelAtmosphere) = [l.tau_5000 for l in atm.layers]

This is a convienince functions for making plots, etc.  Note that it doesn't access quantities in a 
memory-efficient order.
"""
get_tau_5000s(atm::ModelAtmosphere) = [l.tau_5000 for l in atm.layers]
"""
   zs(atm::ModelAtmosphere) = [l.z for l in atm.layers]

This is a convienince functions for making plots, etc.  Note that it doesn't access quantities in a 
memory-efficient order.
"""
get_zs(atm::ModelAtmosphere) = [l.z for l in atm.layers]
"""
   temps(atm::ModelAtmosphere) = [l.temp for l in atm.layers]

This is a convienince functions for making plots, etc.  Note that it doesn't access quantities in a 
memory-efficient order.
"""
get_temps(atm::ModelAtmosphere) = [l.temp for l in atm.layers]
"""
   electron_number_densities(atm::ModelAtmosphere) = [l.electron_number_density for l in atm.layers]

This is a convienince functions for making plots, etc.  Note that it doesn't access quantities in a 
memory-efficient order.
"""
get_electron_number_densities(atm::ModelAtmosphere) = [l.electron_number_density for l in atm.layers]
"""
    number_densities(atm::ModelAtmosphere) = [l.number_density for l in atm.layers]

This is a convienince functions for making plots, etc.  Note that it doesn't access quantities in a 
memory-efficient order.
"""
get_number_densities(atm::ModelAtmosphere) = [l.number_density for l in atm.layers]

"""
    read_model_atmosphere(filename)

Parse the provided model atmosphere file in MARCS ".mod" format.  Returns either a 
`PlanarAtmosphere` or a `ShellAtmosphere`.
"""
function read_model_atmosphere(fname::AbstractString) :: ModelAtmosphere
    open(fname) do f
        #these files are small, so it's not a big deal to load them entirely into memory
        lines = collect(eachline(f)) 
        
        Rind = findfirst(occursin.("adius", lines)) # {rR}adius has uncertain capitalization
        if isnothing(Rind)
            throw(ArgumentError("Can't parse .mod file:  can't detect radius." * 
                                " (should be 1.0 for plane-parallel atmospheres.)"))
        end
        R = parse(Float64, split(lines[Rind])[1])
        planar = R == 1

        i = findfirst(occursin.("Number of depth points", lines))
        if isnothing(i)
            throw(ArgumentError("Can't parse .mod file: can't detect number of layers."))
        end
        nlayers = parse(Int, split(lines[i])[1])

        header = findfirst(occursin.("lgTauR", lines))
        if isnothing(header)
            throw(ArgumentError("Can't parse .mod file: can't find header."))
        end

        layers = map(lines[header+1:header+nlayers]) do line 
            logτ5 = parse(Float64, line[11:17])
            depth = parse(Float64, line[19:28])
            temp = parse(Float64, line[30:36])
            Pe = parse(Float64, line[39:48])
            Pg = parse(Float64, line[49:60])

            # round negative pressures to 0
            Pe = Pe * (Pe > 0)
            Pg = Pg * (Pg > 0)

            nₑ = Pe / (temp*kboltz_cgs) # electron number density
            n = Pg/ (temp*kboltz_cgs)   # total number density

            if planar
                PlanarAtmosphereLayer(10^logτ5, -depth, temp, nₑ, n)
            else
                ShellAtmosphereLayer(10^logτ5, -depth, temp, nₑ, n)
            end
        end

        if planar
            PlanarAtmosphere(layers)
        else
            ShellAtmosphere(layers, R)
        end
     end
end

const _sdss_marcs_atmospheres = let
    path = joinpath(artifact"SDSS_MARCS_atmospheres", "SDSS_MARCS_atmospheres.h5")
    exists = h5read(path, "exists")
    grid = h5read(path, "grid")
    nodes = [h5read(path, "grid_values/$i") for i in 1:5]
    (nodes, exists, grid)
end

"""
    interpolate_marcs(Teff, logg, Fe_M=0, alpha_M=0, C_M=0; kwargs...)
    interpolate_marcs(Teff, logg, A_X; kwargs...)

Returns a model atmosphere obtained by interpolating the MARCS SDSS atmosphere grid, which provides 
atmospheres for varying values of ``T_\\mathrm{eff}``, ``\\log g``, [metals/H], [alpha/metals], 
and [C/metals].  If the `A_X` (a vector of abundances in the format returned by [`format_A_X`](@ref)
and accepted by [`synthesize`](@ref)) is provided instead of `M_H`, `alpha_M`, and `C_M`, the 
solar-relative ratios will be reconstructed assuming Grevesse+ 2007 solar abundances.

The model atmosphere grid is a repacked version of the 
[MARCS SDSS grid](https://dr17.sdss.org/sas/dr17/apogee/spectro/speclib/atmos/marcs/MARCS_v3_2016/Readme_MARCS_v3_2016.txt).

# keyword arguments
- `spherical`: whether or not to return a ShellAtmosphere (as opposed to a PlanarAtmosphere).  By 
  default true when `logg` < 3.5.
- `archive`: The atmosphere archive to use.  This is used to override the default grid for testing.
- `solar_abundances`: (default: `grevesse_2007_solar_abundances`) The solar abundances to use when 
  `A_X` is provided instead of `M_H`, `alpha_M`, and `C_M`. The default is chosen to match that of 
  the atmosphere grid, and is probably no good reason to change it.
- `clamp_abundances`: (default: `false`) allowed when specifying `A_X` direction. Whether or not to 
   clamp the abundance paramerters to be within the range of the MARCS grid to avoid throwing an out 
   of bounds error. Use with caution.
- `perturb_at_grid_values`: when true this will add or a subtract a very small number to each 
   parameter which is exactly at a grid value. This prevents null derivatives, which can cause 
   problems for minimizers.  

!!! warning
    Atmosphere interpolation contributes non-negligeble error to synthesized spectra below 
    Teff ≈ 4250 K. We do not endorse using it for science in that regime. See 
    https://github.com/ajwheeler/Korg.jl/issues/164 for a discussion of the issue.
"""
function interpolate_marcs(Teff, logg, A_X::AbstractVector{<:Real}; 
                           solar_abundances=grevesse_2007_solar_abundances, 
                           clamp_abundances=false, kwargs...)
    M_H = get_metals_H(A_X; solar_abundances=solar_abundances)
    alpha_H = get_alpha_H(A_X; solar_abundances=solar_abundances)
    alpha_M = alpha_H - M_H
    C_H = A_X[6] - solar_abundances[6]
    C_M = C_H - M_H
    if clamp_abundances
        nodes = _sdss_marcs_atmospheres[1]
        M_H = clamp(M_H, nodes[3][1], nodes[3][end])
        alpha_M = clamp(C_M, nodes[4][1], nodes[4][end])
        C_M = clamp(C_M, nodes[5][1], nodes[5][end])
    end
    interpolate_marcs(Teff, logg, M_H, alpha_M, C_M; kwargs...)
end
function interpolate_marcs(Teff, logg, M_H=0, alpha_M=0, C_M=0; spherical=logg < 3.5, 
                           perturb_at_grid_values=false)
    nodes, _, grid = _sdss_marcs_atmospheres
    params = [Teff, logg, M_H, alpha_M, C_M]
    param_names = ["Teff", "log(g)", "[M/H]", "[alpha/M]", "[C/metals]"]
    atm_quants = lazy_multilinear_interpolation(params, nodes, grid, param_names=param_names, 
                                               perturb_at_grid_values=perturb_at_grid_values)

    if spherical
        R = sqrt(G_cgs * solar_mass_cgs / 10^(logg)) 
        ShellAtmosphere(ShellAtmosphereLayer.(atm_quants[:, 4], 
                                              sinh.(atm_quants[:, 5]), 
                                              atm_quants[:, 1],
                                              exp.(atm_quants[:, 2]), 
                                              exp.(atm_quants[:, 3])), R)
    else
        PlanarAtmosphere(PlanarAtmosphereLayer.(atm_quants[:, 4], 
                                               sinh.(atm_quants[:, 5]), 
                                               atm_quants[:, 1],
                                               exp.(atm_quants[:, 2]), 
                                               exp.(atm_quants[:, 3])))
    end
end