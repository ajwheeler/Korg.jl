"""
    read_model_atmosphere(fname)

Parse the provided model atmosphere file in 
["Kurucz" format](https://marcs.astro.uu.se/krz_format.html).
"""
function read_model_atmosphere(fname::AbstractString)
    open(fname) do f
        #these files are small, so it's not a big deal to load them entirely into memory
        lines = collect(eachline(f)) 
        map(lines[14:end]) do line
            toks = strip.(split(line, ','))
            #spherical atmospheres have a sixth column with distance from reference height.
            #ignore that for now.
            (tau=parse(Float64, toks[1]),              #optical depth
             temp=parse(Float64, toks[2]),             #K
             electron_density=parse(Float64, toks[3]), #cm^-3
             number_density=parse(Float64, toks[4]),   #of everything except elections, cm^-3
             density=parse(Float64, toks[5]))          #g cm^-3
        end
    end
end
