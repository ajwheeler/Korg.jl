abstract type ModelAtmosphere end

struct PlanarAtmosphereLayer{F}
    colmass::F                 #g
    temp::F                    #K
    electron_number_density::F #cm^-3
    number_density::F          #cm^-3
    density::F                 #g cm^-3
end

struct PlanarAtmosphere{F} <: ModelAtmosphere
    layers::Vector{PlanarAtmosphereLayer{F}}
end

struct ShellAtmosphereLayer{F}
    colmass::F                 #g
    temp::F                    #K
    electron_number_density::F #cm^-3
    number_density::F          #cm^-3
    density::F                 #g cm^-3
    r::F                       #cm
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
