using Korg: RydbergH_eV

function parse_NORAD_H_cross_sections(filename, norad_format=false)
    cross_sections = []
    lines = readlines(filename)

    #points to "current" line. This could be trivially rewritten to stream line-by-line from disk.
    i = 1

    #skip the big header in norad formatted files
    while !startswith(lines[i], "---------------------------------------------------------")
        i += 1
    end
    i += 1

    i += 1 # the first line tells you which species the file corresponds to, skip to the second

    while i <= length(lines)
        # parse the header specifying the state to get the angular momentum and primary quantum nums
        L = parse(UInt8, lines[i][9:10])
        n = parse(UInt8, lines[i][19:20])

        # the next line tells you the binding energy (in Rydberg energies) and how many energies the 
        # cross-section is evaluated at
        toks = split(lines[i+1])
        binding_E = parse(Float64, toks[1])
        npoints = parse(Int, toks[2])

        # map the electron state the cross-section interpolator in the cross_sections dict.
        Es = Vector{Float64}(undef, npoints)
        σs = Vector{Float64}(undef, npoints)

        i += 2 # Make i point at the first real line of the cross-sections vals
        for j in 1:npoints
            Es[j] = parse(Float64, lines[i+j-1][3:14])
            σs[j] = parse(Float64, lines[i+j-1][26:36])
        end

        push!(cross_sections, (n, L, binding_E * RydbergH_eV, Es * RydbergH_eV, σs))

        i += npoints #move i to the next electron state
    end

    cross_sections
end
