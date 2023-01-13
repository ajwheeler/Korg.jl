# used when downloading model atmosphere archive
using ProgressMeter: Progress, update!, finish!
using Downloads: download
 
abstract type ModelAtmosphere end

struct PlanarAtmosphereLayer{F}
    tau_5000::F                #dimensionless (used for legacy radiative transfer)
    z::F                       #cm
    temp::F                    #K
    electron_number_density::F #cm^-3
    number_density::F          #cm^-3
end

struct PlanarAtmosphere{F} <: ModelAtmosphere
    layers::Vector{PlanarAtmosphereLayer{F}}
end

struct ShellAtmosphereLayer{F}
    tau_5000::F                #dimensionless (used for legacy radiative transfer)
    z::F                       #cm
    temp::F                    #K
    electron_number_density::F #cm^-3
    number_density::F          #cm^-3
end

struct ShellAtmosphere{F} <: ModelAtmosphere
    layers::Vector{ShellAtmosphereLayer{F}}
    R::F #the radius of the star where τ_ros == 1, i.e. the photosphere (not the top)
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

!!! warning
    Korg does not yet support cool (``\\lesssim`` 3500 K) stars.  While it will happily parse their
    model atmospheres, it is very likely to crash when you feed them into `synthesize`.
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
            n = Pg/ (temp*kboltz_cgs)   # non-electron number density

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

# this isn't a const because the model atmosphere doesn't get loaded into memory until 
# interpolate_marcs is called for the first time
global _atmosphere_archive = nothing 

"""
Returns the local where Korg's large (too big for git) data files are stored.  At present, this is
only the model atmosphere archive used by [`interpolate_marcs`](@ref).
"""
function _korg_data_dir()
    if "KORG_DATA_DIR" in keys(ENV)
        joinpath(ENV["KORG_DATA_DIR"])
    else
        joinpath(homedir(), ".korg")
    end
end

"""
Returns the model atmosphere grid, first loading it into memory if necessary. It's not loaded 
at import time because it's large, and model atmosphere interpolation isn't always needed. This
is called automatically when running `interpolate_marcs`.
"""
function get_atmosphere_archive()
    if !isnothing(_atmosphere_archive)
        return _atmosphere_archive
    end
    @info "loading the model atmosphere grid into memory. This will take a few seconds, but will only happen once per julia session."
    path = joinpath(_korg_data_dir(), "SDSS_MARCS_atmospheres.h5")
    if !isfile(path)
        msg = """
        Could not find the model atmosphere archive.  If this is the first time you are running 
        interpolate_marcs, you will need to download the archive by running 

            Korg.download_atmosphere_archive()

        This will download the ~370 MB archive file and place it in ~/.korg/  If you would like 
        to store it somewhere else, you can specify a location with the KORG_DATA_DIR environment 
        variable.
        """
        throw(ErrorException(msg))
    end

    exists = h5read(path, "exists")
    grid = h5read(path, "grid")
    nodes = [h5read(path, "grid_values/$i") for i in 1:5]
    global _atmosphere_archive = (nodes, exists, grid)
end

"""
    download_atmosphere_archive()

Download the data used by [`interpolate_marcs`](@ref), a repacked version of the 
[MARCS SDSS grid](https://dr17.sdss.org/sas/dr17/apogee/spectro/speclib/atmos/marcs/MARCS_v3_2016/Readme_MARCS_v3_2016.txt).
By default, the archive is stored at `.korg/SDSS_MARCS_atmospheres.h5`.  This location can be set 
with the `KORG_DATA_DIR` environment variable.
"""
function download_atmosphere_archive(url="https://korg-data.s3.amazonaws.com/SDSS_MARCS_atmospheres.h5")
    prog = Progress(100, output=stdout, desc="Downloading model atmosphere archive: ")
    finished = false
    function update_progress(total, now)
        # deal with initial total and downloaded sizes being indeterminate
        if now == 0 || isinf(now) || isinf(total) 
            now, total = 0,1 
        end
        # Downloads will call this function a few times with now==total, but we only want to
        # update the progress bar once, to avoid duplicates bars.
        if finished == false
            update!(prog, Int(floor(now/total * 100)))
        end
        if now == total
            finished = true
        end
    end
    data_dir = _korg_data_dir()
    if !isdir(data_dir)
        @info "creating $data_dir"
        mkdir(data_dir)
    end
    download(url, joinpath(data_dir, "SDSS_MARCS_atmospheres.h5"), progress=update_progress);
end

"""
    interpolate_marcs(archive, Teff, logg, Fe_H=0, alpha_Fe=0, C_Fe=0)
    interpolate_marcs(archive, Teff, logg, A_X)

Returns a model atmosphere obtained by interpolating the atmosphere grid `archive`. 
If the `A_X` (a vector of abundances in the format returned by [`format_A_X`](@ref) and accepted by 
[`synthesize`](@ref)) is provided instead of `Fe_H`, `alpha_Fe`, and `C_Fe`, the solar-relative 
ratios will be reconstructed assuming Grevesse+ 2007 solar abundances, with Mg determining the alpha
ratio.

The model atmosphere grid is a repacked version of the 
[MARCS SDSS grid](https://dr17.sdss.org/sas/dr17/apogee/spectro/speclib/atmos/marcs/MARCS_v3_2016/Readme_MARCS_v3_2016.txt).
Before you use `interpolate_marcs` for the first time, you will have to download the grid with by 
running `download_atmosphere_archive()`.
By default, the archive is stored at `.korg/SDSS_MARCS_atmospheres.h5`.  This location can be set 
with the `KORG_DATA_DIR` environment variable.

# keyword arguments
- `spherical`: whether or not to return a ShellAtmosphere (as opposed to a PlanarAtmosphere).  By 
  default true when `logg` < 3.5.
- `archive`: The atmosphere archive to use. For testing purposes.

!!! warning
    Atmosphere interpolation is in beta.  Use with caution as it may be inaccurate.

"""
function interpolate_marcs(Teff, logg, A_X::Vector; kwargs...)
    Fe_H = A_X[26] - grevesse_2007_solar_abundances[26], A[X]
    alpha_Fe = A_X[12]/A_X[26] - grevesse_2007_solar_abundances[12]/grevesse_2007_solar_abundances[26]
    C_Fe = A_X[6]/A_x[26] - grevesse_2007_solar_abundances[6]/grevesse_2007_solar_abundances[26]
    interpolate_marcs(Teff, logg, Fe_H, alpha_Fe, C_Fe; kwargs...)
end
function interpolate_marcs(Teff, logg, Fe_H=0, alpha_Fe=0, C_Fe=0; spherical=logg < 3.5, 
                           archive=get_atmosphere_archive())
    nodes, exists, grid = archive

    params = [Teff, logg, Fe_H, alpha_Fe, C_Fe]
    param_names = ["Teff", "log(g)", "[Fe/H]", "[alpha/Fe]", "[C/Fe]"]
    
    upper_vertex = map(zip(params, param_names, nodes)) do (p, p_name, p_nodes)
        if !(p_nodes[1] <= p <= p_nodes[end])
            throw(ArgumentError("Can't interpolate atmosphere grid.  $(p_name) is out of bounds. ($(p) ∉ [$(first(p_nodes)), $(last(p_nodes))])"))
        end
        findfirst(p .<= p_nodes)
    end
    isexact = params .== getindex.(nodes, upper_vertex) #which params are on grid points?
    
    # allocate 2^n cube for each quantity
    dims = Tuple(2 for _ in upper_vertex) #dimensions of 2^n hypercube
    structure_type = typeof(promote(Teff, logg, Fe_H, alpha_Fe, C_Fe)[1])
    structure = Array{structure_type}(undef, (56, 5, dims...))
     
    #put bounding atmospheres in 2^n cube
    for I in CartesianIndices(dims)
        local_inds = collect(Tuple(I))
        atm_inds = copy(local_inds)
        atm_inds[isexact] .= 2 #use the "upper bound" as "lower bound" if the param is on a grid point
        atm_inds .+= upper_vertex .- 2

        if !exists[atm_inds...] #return nothing if any required nodes don't exist
            @info str(getindex.(nodes, atm_inds), " doesn't exist") 
            return 
        end
        
        structure[:, :, local_inds...] .= grid[:, :, atm_inds...]
    end

    for i in 1:length(params) #loop over Teff, logg, etc.
        isexact[i] && continue #no need to do anything for exact params
        
        # the bounding values of the parameter you are currently interpolating over 
        p1 = nodes[i][upper_vertex[i]-1]
        p2 = nodes[i][upper_vertex[i]]
        
        # inds1 and inds2 are the expressions for the slices through the as-of-yet 
        # uninterpolated quantities (temp, logPg, etc) for each node value of the
        # quantity being interpolated
        # inds1 = (1, 1, 1, ..., 1, 1, :, :, ...)
        # inds2 = (1, 1, 1, ..., 1, 2, :, :, ...)
        inds1 = vcat([1 for _ in 1:i-1], 1, [Colon() for _ in i+1:length(params)])
        inds2 = vcat([1 for _ in 1:i-1], 2, [Colon() for _ in i+1:length(params)])
        
        x = (params[i] - p1) / (p2 - p1) #maybe try using masseron alpha later
        for structure_ind in 1:5 #TODO fix
            structure[:, structure_ind, inds1...] = (1-x)*structure[:, structure_ind, inds1...] + x*structure[:, structure_ind, inds2...]
        end
    end
   
    atm_quants = (structure[:, :, ones(Int, length(params))...])
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