abstract type ModelAtmosphere end

struct PlanarAtmosphereLayer{F}
    z::F                       #g
    temp::F                    #K
    electron_number_density::F #cm^-3
    number_density::F          #cm^-3
    density::F                 #g cm^-3
end

struct PlanarAtmosphere{F} <: ModelAtmosphere
    layers::Vector{PlanarAtmosphereLayer{F}}
end

struct ShellAtmosphereLayer{F}
    r::F                       #cm
    temp::F                    #K
    electron_number_density::F #cm^-3
    number_density::F          #cm^-3
    density::F                 #g cm^-3
end

struct ShellAtmosphere{F} <: ModelAtmosphere
    layers::Vector{ShellAtmosphereLayer{F}}
end

"""
    PlanarAtmosphere(atm::ShellAtmosphere)

Construct a planar atmosphere with the data from a shell atmosphere.  Mostly useful for testing.
"""
function PlanarAtmosphere(atm::ShellAtmosphere)
    PlanarAtmosphere([PlanarAtmosphereLayer(l.colmass, l.temp, l.electron_number_density, 
                                            l.number_density, l.density) for l in atm.layers])
end

"""
    ShellAtmosphere(atm::PlanarAtmosphere, R)

Construct a shell atmosphere with the data from a planar atmosphere and an outer radius.  Mostly 
useful for testing.
"""
function ShellAtmosphere(atm::PlanarAtmosphere, R)
    Δcolmass = diff((l->l.colmass).(atm.layers))
    Δs = 0.5([0 ; Δcolmass] + [Δcolmass; Δcolmass[end]]) ./ (l->l.density).(atm.layers)
    rs = R .- cumsum(Δs)

    ShellAtmosphere([ShellAtmosphereLayer(l.colmass, l.temp, l.electron_number_density, 
                                          l.number_density, l.density, r) 
                     for (l, r) in zip(atm.layers, rs)])
end

#pretty-printing
function Base.show(io::IO, m::MIME"text/plain", atm::A) where A <: ModelAtmosphere
    print(io, "$(A) with $(length(atm.layers)) layers")
end

"""
    read_model_atmosphere(filename; truncate_at_10000K=true)

Parse the provided model atmosphere file in 
["Kurucz" format](https://marcs.astro.uu.se/krz_format.html). Returns either a `PlanarAtmosphere` or
a `ShellAtmosphere`.

When `truncate_at_10000K` is true, layers with temperatures greater than 10,000 Kelvin will be 
elided.  These are typically at the deepest level, where optical depth is very large, and have 
minimal impact on the surface spectrum.  They are ignored by default because Korg's default 
partition and molecular equillibrium functions are tabulated only up to that tempurature.
"""
function read_model_atmosphere(fname::AbstractString; truncate_at_10000K=true)
    open(fname) do f
        #these files are small, so it's not a big deal to load them entirely into memory
        lines = collect(eachline(f)) 

        toks = strip.(split(lines[2], ' '))
        spherical = "RADIUS=" in toks
        R_header = spherical ? parse(Float64, toks[end-1]) : NaN

        layers = map(lines[14:end]) do line
            toks = strip.(split(line, ','))
            if toks[end] == ""
                toks = toks[1:end-1]
            end
            numbers = parse.(Float64, toks)
            if spherical
                #spherical shell atmospheres have a sixth column with distance from reference height
                numbers[6] += R_header
                ShellAtmosphereLayer(numbers...)
            else
                PlanarAtmosphereLayer(numbers...)
            end
         end

        if truncate_at_10000K
            filter!(layers) do layer
                layer.temp < 10_000
            end
        end

         if spherical
             ShellAtmosphere(layers)
         else
             PlanarAtmosphere(layers)
         end
    end
end

function read_mod_atmosphere(fname::AbstractString; truncate_at_10000K=true) :: ModelAtmosphere
    open(fname) do f
        #these files are small, so it's not a big deal to load them entirely into memory
        lines = collect(eachline(f)) 
        
        Rind = findfirst(occursin.("Radius", lines))
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
            temp = nums[5]
            Pe = nums[6]
            Pg = nums[7]

            if planar
                PlanarAtmosphereLayer(-nums[4], temp, Pe/(temp*kboltz_cgs), Pg/(temp*kboltz_cgs), 
                                      nums[13])
            else
                ShellAtmosphereLayer(R-nums[4], temp, Pe/(temp*kboltz_cgs), Pg/(temp*kboltz_cgs), 
                                      nums[13])
            end
        end

        if truncate_at_10000K
            filter!(layers) do layer
                layer.temp < 10_000
            end
        end

        #TODO handle spherical atmosphers
        if planar
            PlanarAtmosphere(layers)
        else
            ShellAtmosphere(layers)
        end
    end
end
