# used when downloading model atmosphere archive
using ProgressMeter: Progress, update!, finish!
using Pkg.Artifacts: @artifact_str
import Interpolations
 
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
    get_tau_5000s(atm::ModelAtmosphere) = [l.tau_5000 for l in atm.layers]

This is a convienince functions for making plots, etc.  Note that it doesn't access quantities in a 
memory-efficient order.
"""
get_tau_5000s(atm::ModelAtmosphere) = [l.tau_5000 for l in atm.layers]
"""
   get_zs(atm::ModelAtmosphere) = [l.z for l in atm.layers]

This is a convienince functions for making plots, etc.  Note that it doesn't access quantities in a 
memory-efficient order.
"""
get_zs(atm::ModelAtmosphere) = [l.z for l in atm.layers]
"""
   get_temps(atm::ModelAtmosphere) = [l.temp for l in atm.layers]

This is a convienince functions for making plots, etc.  Note that it doesn't access quantities in a 
memory-efficient order.
"""
get_temps(atm::ModelAtmosphere) = [l.temp for l in atm.layers]
"""
   get_electron_number_densities(atm::ModelAtmosphere) = [l.electron_number_density for l in atm.layers]

This is a convienince functions for making plots, etc.  Note that it doesn't access quantities in a 
memory-efficient order.
"""
get_electron_number_densities(atm::ModelAtmosphere) = [l.electron_number_density for l in atm.layers]
"""
    get_number_densities(atm::ModelAtmosphere) = [l.number_density for l in atm.layers]

This is a convienince functions for making plots, etc.  Note that it doesn't access quantities in a 
memory-efficient order.
"""
get_number_densities(atm::ModelAtmosphere) = [l.number_density for l in atm.layers]
"""
    get_gas_pressures(atm::ModelAtmosphere) 

This is a convienince functions for making plots, etc.  Note that it doesn't access quantities in a
memory-efficient order.
"""
get_gas_pressures(atm) = [l.number_density * kboltz_cgs * l.temp for l in atm.layers]

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

# used for the standard and low-metallicity grids.  Lazy linear interp is used for these, so it's 
# just a matter of reading the data and setting up the mmap.
function _prepare_linear_atmosphere_archive(path)
    h5open(path, "r") do f
        grid = HDF5.readmmap(f["grid"])
        params = read(f["grid_parameter_names"])
        @assert params == ["Teff", "logg", "metallicity", "alpha", "carbon"][1:length(params)]
        nodes = [read(f["grid_values/$i"]) for i in 1:length(params)]
        nodes, grid
    end
end
_sdss_marcs_atmospheres = let 
    # note to self: don't put files in a directory before you tarball it next time.  It's redundant!
    path = joinpath(artifact"SDSS_MARCS_atmospheres_v2", "SDSS_MARCS_atmospheres", "SDSS_MARCS_atmospheres.h5")
    _prepare_linear_atmosphere_archive(path)
end
_low_Z_marcs_atmospheres = let 
    path = joinpath(artifact"MARCS_metal_poor_atmospheres", "MARCS_metal_poor_atmospheres", "MARCS_metal_poor_atmospheres.h5")
    _prepare_linear_atmosphere_archive(path)
end

# cubic interp is used for the cool dwarfs, so we need to set up the interpolator. This takes more 
# cpu/memory.
function _prepare_cool_dwarf_atm_archive(grid, nodes)
    nodes_ranges = [range(first(n), last(n), length(n)) for n in nodes]   
    @assert all(nodes_ranges .== nodes)

    nlayers = size(grid, 1)
    knots = tuple(1f0:nlayers, 1f0:5f0, nodes_ranges...)

    itp = Interpolations.scale(Interpolations.interpolate(
        grid, (Interpolations.NoInterp(),
               Interpolations.NoInterp(),
               Interpolations.BSpline(Interpolations.Cubic()),
               Interpolations.BSpline(Interpolations.Cubic()), 
               Interpolations.BSpline(Interpolations.Cubic()),
               Interpolations.BSpline(Interpolations.Cubic()),
               Interpolations.BSpline(Interpolations.Cubic()))),
          knots)
    itp, nlayers
end
const _cool_dwarfs_atm_itp = let 
    path = joinpath(artifact"resampled_cool_dwarf_atmospheres", "resampled_cool_dwarf_atmospheres", "resampled_cool_dwarf_atmospheres.h5")
    grid, nodes = h5open(path, "r") do f
        read(f["grid"]), [read(f["grid_values/$i"]) for i in 1:5]
    end
    _prepare_cool_dwarf_atm_archive(grid, nodes)
end

"""
    interpolate_marcs(Teff, logg, A_X; kwargs...)
    interpolate_marcs(Teff, logg, m_H=0, alpha_m=0, C_m=0; kwargs...)

Returns a model atmosphere computed by interpolating models from [MARCS](https://marcs.astro.uu.se/)
((Gustafsson+ 2008)[https://ui.adsabs.harvard.edu/abs/2008A&A...486..951G/abstract]).
Along with `Teff` and `logg`, the atmosphere is specified by `m_H`, `alpha_m`, and `C_m`, which can 
be automatically determined from an `A_X` abundance vector (the recommended method, 
see [`format_A_X`](@ref)). Note that the MARCS atmosphere models were constructed with the 
Grevesse+ 2007 solar abundances (`Korg.grevesse_2007_solar_abundances`). This is handled 
automatically when `A_X` is provided.

`interpolate_marcs` uses three different interpolation schemes for different stellar parameter
regimes. In the standard case the model atmosphere grid is [the one generated for
SDSS](https://dr17.sdss.org/sas/dr17/apogee/spectro/speclib/atmos/marcs/MARCS_v3_2016/Readme_MARCS_v3_2016.txt),
transformed and linearly interpolated. For cool dwarfs (`Teff` ≤ 4000 K, `logg` ≥ 3.5), the grid is
resampled onto unchanging `tau_5000` values and interpolated with a cubic spline. For
low-metallicity models (-5 ≤ `m_H` < -2.5), a grid of standard composition (i.e. fixed alpha and C)
atmospheres is used.  (The microturbulence is 1km/s for dwarfs and 2km/s for giants and the mass for
spherical models is 1 solar mass.) The interpolation method is the same as in the standard case. See
[Wheeler+ 2024](https://ui.adsabs.harvard.edu/abs/2023arXiv231019823W/abstract) for more details and
a discussion of errors introduced by model atmosphere interpolation. (Note that the cubic scheme for
cool dwarfs is referred to as not-yet-implemented in the paper but is now available.)

# keyword arguments
- `spherical`: whether or not to return a ShellAtmosphere (as opposed to a PlanarAtmosphere).  By 
  default true when `logg` < 3.5.
- `solar_abundances`: (default: `grevesse_2007_solar_abundances`) The solar abundances to use when 
  `A_X` is provided instead of `M_H`, `alpha_M`, and `C_M`. The default is chosen to match that of 
  the atmosphere grid, and if you change it you are likely trying to do something else.
- `clamp_abundances`: (default: `false`) allowed only when specifying `A_X`. Whether or not to 
   clamp the abundance parameters to be within range to avoid throwing an out of bounds error.
- `perturb_at_grid_values` (default: `true`): whether or not to add or a subtract a very small 
   number to each parameter which is exactly at a grid value. This prevents null derivatives, which 
   can cause problems for minimizers. 
- `resampled_cubic_for_cool_dwarfs` (default: `true`): whether or not to used specialized method for cool dwarfs.
- `archives`: A tuple containing the atmosphere grids to use.  For testing purposes.
"""
function interpolate_marcs(Teff, logg, A_X::AbstractVector{<:Real}; 
                           solar_abundances=grevesse_2007_solar_abundances, 
                           clamp_abundances=false, 
                           archives=(_sdss_marcs_atmospheres, _cool_dwarfs_atm_itp, _low_Z_marcs_atmospheres),
                           kwargs...)

    m_H = get_metals_H(A_X; solar_abundances=solar_abundances, 
                       alpha_elements=[6 ; default_alpha_elements]) # ignore C in addition to alphas
    alpha_H = get_alpha_H(A_X; solar_abundances=solar_abundances)
    alpha_m = alpha_H - m_H
    C_H = A_X[6] - solar_abundances[6]
    C_m = C_H - m_H

    if clamp_abundances
        standard_nodes = archives[1][1]
        low_Z_nodes = archives[3][1]
        m_H = clamp(m_H, low_Z_nodes[3][1], standard_nodes[3][end])
        alpha_m = clamp(alpha_m, standard_nodes[4][1], standard_nodes[4][end])
        C_m = clamp(C_m, standard_nodes[5][1], standard_nodes[5][end])
    end

    if m_H < -2.5
        # these are the only values allowed for low-metallicity models
        alpha_m = 0.4
        C_m = 0
    end

    interpolate_marcs(Teff, logg, m_H, alpha_m, C_m; archives=archives, kwargs...)
end
function interpolate_marcs(Teff, logg, m_H=0, alpha_m=0, C_m=0; spherical=logg < 3.5, 
                           perturb_at_grid_values=true, resampled_cubic_for_cool_dwarfs=true,
                           archives=(_sdss_marcs_atmospheres, _cool_dwarfs_atm_itp, _low_Z_marcs_atmospheres))
    # cool dwarfs
    if Teff <= 4000 && logg >= 3.5 && m_H >= -2.5 && resampled_cubic_for_cool_dwarfs
        itp, nlayers = archives[2]
        atm_quants = itp(1:nlayers, 1:5, Teff, logg, m_H, alpha_m, C_m)
        PlanarAtmosphere(PlanarAtmosphereLayer.(
            atm_quants[:, 4],
            sinh.(atm_quants[:, 5]), 
            atm_quants[:, 1],
            exp.(atm_quants[:, 2]), 
            exp.(atm_quants[:, 3])))
    else
        # low metallicity
        if m_H < -2.5
            nodes, grid = archives[3]
            if alpha_m != 0.4 || C_m != 0
                throw(ArgumentError("For low metallicities ([m_H < -2.5]), it is required that alpha_M = 0.4 and C_M = 0"))
            end
            params = [Teff, logg, m_H]
            param_names = ["Teff", "log(g)", "[M/H]"]
        # standard 
        else
            nodes, grid = archives[1]
            params = [Teff, logg, m_H, alpha_m, C_m]
            param_names = ["Teff", "log(g)", "[M/H]", "[alpha/M]", "[C/metals]"]
        end

        atm_quants = lazy_multilinear_interpolation(params, nodes, grid, param_names=param_names, 
                                                    perturb_at_grid_values=perturb_at_grid_values)

        # grid atmospheres are allowed to have to NaNs to represent layers that should be dropped. 
        nanmask = .! isnan.(atm_quants[:, 4]) # any column will do. This is τ_5000.

        if spherical
            R = sqrt(G_cgs * solar_mass_cgs / 10^(logg)) 
            ShellAtmosphere(ShellAtmosphereLayer.(atm_quants[nanmask, 4], 
                                                  sinh.(atm_quants[nanmask, 5]), 
                                                  atm_quants[nanmask, 1],
                                                  exp.(atm_quants[nanmask, 2]), 
                                                  exp.(atm_quants[nanmask, 3])), R)
        else
            PlanarAtmosphere(PlanarAtmosphereLayer.(atm_quants[nanmask, 4], 
                                                   sinh.(atm_quants[nanmask, 5]), 
                                                   atm_quants[nanmask, 1],
                                                   exp.(atm_quants[nanmask, 2]), 
                                                   exp.(atm_quants[nanmask, 3])))
        end
    end
end
# handle the case where Teff, logg, and [m/H] are integers. As long as not all (interpolated) params 
# are passed in as integers, there's no problem. 
interpolate_marcs(Teff::Int, logg::Int, m_H::Int, args...; kwargs...) = interpolate_marcs(Float64(Teff), args...; kwargs...)