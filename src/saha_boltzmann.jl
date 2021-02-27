function setup_ionization_energies(fname=joinpath(_data_dir, 
                                                  "BarklemCollet2016-ionization_energies.dat"))
    open(fname, "r") do f
        d = Dict{String, Vector{Float64}}()
        for line in eachline(f)
            if line[1] != '#'        
                toks = split(strip(line))
                #the first token is the atomiz number, which we ignore
                d[toks[2]] = parse.(Float64, toks[3:end])
            end
        end
        d
    end
end

function saha(species, T, upto=3)
     
end
