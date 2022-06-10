abstract type ModelAtmosphere end

struct PlanarAtmosphereLayer{F}
    tau_5000::F                #dimensionless
    z::F                       #g
    temp::F                    #K
    electron_number_density::F #cm^-3
    number_density::F          #cm^-3
end

struct PlanarAtmosphere{F} <: ModelAtmosphere
    layers::Vector{PlanarAtmosphereLayer{F}}
end

struct ShellAtmosphereLayer{F}
    tau_5000::F                #dimensionless
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
        blockA = lines[header+1:header+nlayers]
        blockB = lines[header+nlayers+2 : header+1+2nlayers]

        layers = map(zip(blockA, blockB)) do (lineA, lineB)
            nums = parse.(Float64, strip.(split(lineA * lineB)))
            temp = nums[5] # five digits doesn't let you go higher than our partition funcs top out.
            nₑ = nums[6]/(temp*kboltz_cgs) #electron number density
            n = nums[7]/(temp*kboltz_cgs)  #non-electron number density
            if planar
                PlanarAtmosphereLayer(10^nums[3], -nums[4], temp, nₑ, n)
            else
                ShellAtmosphereLayer(10^nums[3], -nums[4], temp, nₑ, n)
            end
        end

        if planar
            PlanarAtmosphere(layers)
        else
            ShellAtmosphere(layers, R)
        end
     end
end
