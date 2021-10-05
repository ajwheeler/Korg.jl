abstract type ModelAtmosphere end

struct PlanarAtmosphere <: ModelAtmosphere
    layers::Vector{NamedTuple}
end

struct ShellAtmosphere <: ModelAtmosphere
    layers::Vector{NamedTuple}
    R #stellar radius
end

"""
    read_model_atmosphere(filename)

Parse the provided model atmosphere file in 
["Kurucz" format](https://marcs.astro.uu.se/krz_format.html).
"""
function read_model_atmosphere(fname::AbstractString)
    open(fname) do f
        #these files are small, so it's not a big deal to load them entirely into memory
        lines = collect(eachline(f)) 

        toks = strip.(split(lines[2], ' '))
        spherical = "RADIUS=" in toks
        R_header = spherical ? parse(Float64, toks[end-1]) : NaN

        #units:  g cm^-2   K      cm^-3                     cm^-3            g cm^-3
        cols = [:colmass, :temp, :electron_number_density, :number_density, :density] 
        layers = map(lines[14:end]) do line
            toks = strip.(split(line, ','))
            if toks[end] == ""
                toks = toks[1:end-1]
            end
            numbers = parse.(Float64, toks)
            if spherical
                #spherical shell atmospheres have a sixth column with distance from reference height
                numbers[6] += R_header
                (; zip([cols; :r], numbers)...)
            else
                (; zip(cols, numbers)...)
            end
         end
         if spherical
             ShellAtmosphere(layers, NaN)
         else
             PlanarAtmosphere(layers)
         end
    end
end
