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

#TODO add tests
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

module InterpolateSDSSMARCS 
using HDF5: h5read
using ..Korg: PlanarAtmosphere, PlanarAtmosphereLayer, ShellAtmosphere, ShellAtmosphereLayer

atmosphere_archive = "../../korg_files/atmospheres/SDSS_MARCS_atmospheres.h5"

#load from file
planar_exists = h5read(atmosphere_archive, "planar/exists")
planar_grid = h5read(atmosphere_archive, "planar/grid")
planar_nodes = [h5read(atmosphere_archive, "planar/grid_values/$i") for i in 1:5]
spherical_exists = h5read(atmosphere_archive, "spherical/exists")
spherical_grid = h5read(atmosphere_archive, "spherical/grid")
R_grid = h5read(atmosphere_archive, "spherical/R_grid")
spherical_nodes = [h5read(atmosphere_archive, "spherical/grid_values/$i") for i in 1:5]

function interpolate_marcs(Teff, logg, metallicity=0, alpha=0, carbon=0;
        spherical=logg < 3.5,
        nodes=(spherical ? spherical_nodes : planar_nodes),
        exists=(spherical ? spherical_exists : planar_exists),
        grid=(spherical ? spherical_grid : planar_grid),
        R_grid=(spherical ? R_grid : nothing))

    #TODO calculate R using M = 1
    
    params = [Teff, logg, metallicity, alpha, carbon]
    
    upper_vertex = map(zip(params, nodes)) do (p, p_nodes)
        @assert p_nodes[1] <= p <= p_nodes[end]
        findfirst(p .<= p_nodes)
    end
    isexact = params .== getindex.(nodes, upper_vertex) #which params are on grid points?
    
    # allocate 2^n cube for each quantity
    dims = Tuple(2 for _ in upper_vertex) #dimensions of 2^n hypercube
    structure_type = typeof(promote(Teff, logg, metallicity, alpha, carbon)[1])
    structure = Array{structure_type}(undef, (56, 5, dims...))
    if spherical
        R = Array{structure_type}(undef, dims)
    else
        R = [0]
    end
     
    #put bounding atmospheres in 2^n cube
    for I in CartesianIndices(dims)
        local_inds = collect(Tuple(I))
        atm_inds = copy(local_inds)
        atm_inds[isexact] .= 2 #use the "upper bound" as "lower bound" if the param is on a grid point
        atm_inds .+= upper_vertex .- 2

        if !exists[atm_inds...] #return nothing if any required nodes don't exist
            println(getindex.(nodes, atm_inds), " doesn't exist") 
            return 
        end
        
        structure[:, :, local_inds...] .= grid[:, :, atm_inds...]
        spherical && (R[local_inds...] = R_grid[atm_inds...])
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
        if spherical
            R[inds1...] = (1-x)*R[inds1...] + x*R[inds2...]
        end
    end
    
    Float64.(structure[:, :, ones(Int, length(params))...]), Float64(R[1])
end

function assemble_atm((m, R))
    m, R = Float64.(m), Float64(R)
    if R == 0
        PlanarAtmosphere(PlanarAtmosphereLayer.(m[:, 4], sinh.(m[:, 5]), m[:, 1],
                                                          exp.(m[:, 2]), exp.(m[:, 3])))
    else
        ShellAtmosphere(ShellAtmosphereLayer.(m[:, 4], sinh.(m[:, 5]), m[:, 1],
                                                        exp.(m[:, 2]), exp.(m[:, 3])), R)
    end
end

end