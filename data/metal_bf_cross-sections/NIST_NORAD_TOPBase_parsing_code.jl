using Base
using Interpolations: LinearInterpolation, deduplicate_knots!
using Korg: partition_funcs, Species, ismolecule, get_atoms, _data_dir, move_bounds, @species_str
using Korg: hplanck_eV, c_cgs, RydbergH_eV, kboltz_eV # constants
using CSV, DataFrames #for NIST energy level parsing

struct ElectronState
    spin_multiplicity::UInt8 # 2S + 1
    L::UInt8                 # total orbital angular momentum 
    P::UInt8                 # parity 0->even, 1->odd
    iLV::UInt8               # index of the level within its spectroscopic series 
    function ElectronState(spin_multiplicity, L, P, iLV)
        @assert 1 <= spin_multiplicity <= 8 #TODO explain
        @assert P == 0 || P == 1            # parity must be either even or odd
        @assert iLV > 0                     # these are 1-indexed
        new(spin_multiplicity, L, P, iLV)
    end
end

function Base.show(io::IO, state::ElectronState)
    orbital_letter = if state.L == 0
        "S"
    elseif state.L == 1
        "P"
    else
        "DFGHIJKLMNO"[state.L - 1] #2->D, 3->F, etc.
    end
    parity_char = if state.P == 0
        ""
    else
        "*"
    end
    print(io, state.spin_multiplicity, orbital_letter, parity_char, " (level ", state.iLV, ")")
end

function statistical_weight(state::ElectronState)
    state.spin_multiplicity * (2*state.L + 1)
end


const _FINAL_LINE =  "    0    0    0    0"

# In a couple tables, there's a place where the photon energy is listed out of order. This seems to
# be okay given that this is in a section of the table where the cross-section is zero.
const _species_with_problematic_subtables = (Species("He_I"), Species("C_I"), Species("C_II"),
                                             Species("Al_II"), Species("O_II"))

function parse_TOPBase_cross_sections(filename, norad_format=false)
    cross_sections = Dict{ElectronState, Any}()
    lines = readlines(filename)

    #points to "current" line. This could be trivially rewritten to stream line-by-line from disk.
    i = 1 

    #skip the big header in norad formatted files
    if norad_format
        while (lines[i] != "-----------------------------------------------------------------------")
            i += 1
        end
        i += 1
    end

    i += 1 # the first line tells you which species the file corresponds to, skip to the second

    while (i <= length(lines)) && (!startswith(lines[i], _FINAL_LINE))
        #parse the header specifying the state
        state = ElectronState(parse(UInt8, lines[i][4:5]), parse(UInt8, lines[i][9:10]),
                              parse(UInt8, lines[i][14:15]), parse(UInt8, lines[i][19:20]))

        #the next line tells you how many energies the cross-section is evaluated at
        # for example: "     479    489" indicates 489 points. 
        # (The "479" is an internal thing. I think it's the number of points calculated with the 
        # R-matrix method directly, as opposed to being interpolated.)
        npoints = parse(Int, lines[i+1][9:14])
        # The line after that indicates the binding energy.
        # for example: "  1.000000E+00    0.0100" indicates that the binding energy is 1 Ryd.
        binding_energy = parse(Float32, lines[i+2][3:14])

        # map the electron state the cross-section interpolator in the cross_sections dict.
        Es = Vector{Float64}(undef, npoints)
        σs = Vector{Float64}(undef, npoints)

        i += 3 # Make i point at the first real line of the cross-sections vals
        for j in 1:npoints
            Es[j] = parse(Float32, lines[i+j-1][3:14])
            σs[j] = parse(Float32, lines[i+j-1][16:24])
        end
        deduplicate_knots!(Es, move_knots=true) # shift repeated E values to the next float 
        #itp = LinearInterpolation(Es, σs, extrapolation_bc=0.0)

        cross_sections[state] = binding_energy * RydbergH_eV, Es * RydbergH_eV, σs

        i += npoints #move i to the next electron state
    end

    cross_sections
end


"""
This returns spin_multiplicity, L, parity, but NOT an ElectronState.  We have to figure out the 
order within each term for that.  It's used for parsing the NIST files.
"""
function electron_state(term)
    L_numbers = Dict(['S', 'P', 'D', 'F', 'G', 'H', 'I', 'J'] .=> 0:7)

    has_prefex = !isdigit(term[1])
    i = if has_prefex #the index of the char with the spin multiplicity
        3
    else
        1
    end
    
    spin_multiplicity = parse(Int, term[i])
    L = L_numbers[term[i+1]]
    
    parity = term[end] == '*'
    
    spin_multiplicity, L, parity
end

"""
This returns a dict of ElectronState => energy pairs.
"""
function parse_NIST_energy_levels(data_dir, spec)
    df = CSV.read(joinpath(data_dir , string(spec)*".csv"), DataFrame)
    rename!(df, "Level (eV)"=>"level")
    select!(df, ["Configuration", "Term", "J", "level"])
    for col in names(df)
        df[!, col] = (x -> String(x[3:end-1])).(df[!, col])
    end
    
    filter!(df) do row # remove uncertain levels and limits, and other weird stuff
        ( (row.Term != "") && isdigit(row.Term[1]) && (!in('[', row.Term)) 
            && (!in('?', row.Term)) && (!in('?', row.level)) )
    end 
    df.level = map(df.level) do l
        parse(Float64, if l[1] == '['
                l[2:end-1]
            else
                l
            end
            )
    end

    #get the numbers which specify the state: 2S+1, L, π
    df.term_nos = electron_state.(df.Term)
    #now calculate iLV, the energy order of the states within each term
    combine(groupby(df,:term_nos), sdf -> sort(sdf,:level), :term_nos => eachindex => :iLV)
    #df.iLV .= 0
    combine(groupby(df, :term_nos)) do sdf
        sdf.iLV .= sortperm(sdf.level)
    end
    
    df.electron_state = ElectronState.(first.(df.term_nos), 
                                        (x->x[2]).(df.term_nos),
                                        last.(df.term_nos),
                                        df.iLV)
    
    #make sure the hacky term parsing worked
    @assert all(startswith.(string.(df.electron_state), df.Term))

    e_levels = Dict{ElectronState, Float64}()
    for row in eachrow(df)
        e_levels[row.electron_state] = row.level
    end
    e_levels
end
