using Base
using Statistics: mean
using Interpolations: linear_interpolation, deduplicate_knots!
using Korg
using Korg: Species, ismolecule, get_atoms, _data_dir, @species_str
using Korg: hplanck_eV, c_cgs, RydbergH_eV, kboltz_eV # constants
using CSV, DataFrames #for NIST energy level parsing

struct ElectronState
    spin_multiplicity::UInt8 # 2S + 1
    L::UInt8                 # total orbital angular momentum 
    P::UInt8                 # parity 0->even, 1->odd
    iLV::UInt8               # index of the level within its spectroscopic series 
    function ElectronState(spin_multiplicity, L, P, iLV)
        @assert 1 <= spin_multiplicity <= 8
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
        string.(collect("DFGHIJKLMNO"))[state.L-1] #2->D, 3->F, etc.
    end
    parity_char = if state.P == 0
        ""
    else
        "*"
    end
    print(io, state.spin_multiplicity, orbital_letter, parity_char, " (level ", state.iLV, ")")
end

const _FINAL_LINE = "    0    0    0    0"

# In a couple tables, there's a place where the photon energy is listed out of order. This seems to
# be okay given that this is in a section of the table where the cross-section is zero.
const _species_with_problematic_subtables = (Species("He_I"), Species("C_I"), Species("C_II"),
                                             Species("Al_II"), Species("O_II"))

function parse_TOPBase_cross_sections(filename, norad_format=false)
    cross_sections = Dict{ElectronState,Any}()
    lines = readlines(filename)

    #points to "current" line. This could be trivially rewritten to stream line-by-line from disk.
    i = 1

    #skip the big header in norad formatted files
    if norad_format
        while (lines[i] !=
               "-----------------------------------------------------------------------")
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
        deduplicate_knots!(Es; move_knots=true) # shift repeated E values to the next float 

        cross_sections[state] = binding_energy * RydbergH_eV, Es * RydbergH_eV, σs

        i += npoints #move i to the next electron state
    end

    cross_sections
end

"""
This returns spin_multiplicity, L, parity, but NOT an ElectronState.  We have to figure out the
order within each term for that.  It's used for parsing the NIST files.
"""
function parse_term(term)
    L_numbers = Dict((c for c in "SPDFGHIJKLMNOQRTUVWXYZ") .=> 0:21)

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
This returns a dict of ElectronState => (energy, g) for a single species.
"""
function parse_NIST_energy_levels(path)
    df = CSV.read(path, DataFrame)
    rename!(df, "Level (eV)" => "level")
    select!(df, ["Configuration", "Term", "J", "level"])

    # strip off the =""
    for col in names(df)
        df[!, col] = (x -> String(x[3:end-1])).(df[!, col])
    end

    #parse J, calculate g
    Jmatch = match.(r"(\d+)(\/2)?", df.J)
    df = df[.!isnothing.(Jmatch), :]
    df.J = map(Jmatch[.!isnothing.(Jmatch)]) do m
        denom = isnothing(m[2]) ? 1 : 2
        parse(Int, m[1]) // denom
    end
    df.g = Int.(2df.J .+ 1)

    for i in 1:size(df, 1) #strip off lower-case letter term prefixes
        if (df.Term[i] != "") && (length(df.Term[i]) > 2) && islowercase(df.Term[i][1]) &&
           (df.Term[i][2] == ' ')
            df.Term[i] = df.Term[i][3:end]
        end
    end

    # remove levels for various reasons
    filter!(df) do row
        ((row.Term != "")       # self-explanatory
         && (row.J != "")         # self-explanatory
         && (!in('o', row.J))     # to catch "1 or 2 or 3" etc.
         && (row.level != " ")    # self-explanatory
         && isdigit(row.Term[1])  # This is some  weird term notation that I can't figure out.
         && (!in('[', row.Term))  # This is some  weird term notation that I can't figure out.
         && (!in('?', row.Term))  # why bother
         && any(isuppercase(c) for c in row.Term) # written with non-JK coupling scheme
         && (!in(',', row.J))     # these just represent multiple levels(?). Annoying to parse.
         && (!in('a', row.level)) # autoionizing
         && !in('+', row.level))  # these have been observed but only relative energies are known
    end

    #parse energy level ignoring []'s and ?'s
    df.level = map(df.level) do l
        if (l[1] == '[') || (l[1] == '(') #remove brackets
            l = l[2:end-1]
        end
        if (l[end] == '?') #strip suffixes
            l = l[1:end-1]
        end
        parse(Float64, l)
    end

    #group and combine by unique configuration/term pair
    df.config_term = tuple.(df.Configuration, df.Term)
    #df.iLV .= 0
    gdf = groupby(df, :config_term)
    for group in gdf
        @assert all(group.Term .== group.Term[1])
    end
    df = combine(gdf, :g => sum => :g,
                 :level => mean => :level,
                 :Term => first => :Term,
                 :Configuration => first => :Configuration)

    #get the numbers which specify the state: 2S+1, L, π
    df.term_nos = parse_term.(df.Term)
    #now calculate iLV, the energy order of the states within each term
    combine(groupby(df, :term_nos)) do sdf
        sdf.iLV .= sortperm(sdf.level)
    end

    df.electron_state = ElectronState.(first.(df.term_nos),
                                       (x -> x[2]).(df.term_nos),
                                       last.(df.term_nos),
                                       df.iLV)

    #make sure the hacky term parsing worked
    #inds = findall(.!startswith.(string.(df.electron_state), df.Term))
    #display([df.electron_state[inds] df.Term[inds]])
    #@assert all(startswith.(string.(df.electron_state), df.Term))

    e_levels = Dict{ElectronState,Tuple{Float64,Int}}()
    for row in eachrow(df)
        e_levels[row.electron_state] = (row.level, row.g)
    end
    e_levels
end

"""
Returns cross section in megabarns summed from all electron configurations for a given species.
λs should be in cm, Ts in K.  data_dir should contain the local path to Korg_data
(https://github.com/ajwheeler/Korg_data).
"""
function single_species_bf_cross_section(spec::Species, λs, Ts, data_dir)
    #this method handles loading the correct files if you just specify the species
    @assert !ismolecule(spec)

    Z = get_atoms(spec.formula)[1]
    n_electrons = Z - spec.charge

    filename = "p" * lpad(Z, 2, "0") * "." * lpad(n_electrons, 2, "0") * ".dat"
    #use NORAD tables when available, TOPBase otherwise
    path, norad_format = if isfile(joinpath(data_dir, "bf_cross_sections", "NORAD", filename))
        joinpath(data_dir, "bf_cross_sections", "NORAD", filename), true
    else
        joinpath(data_dir, "bf_cross_sections", "TOPBase", filename), false
    end
    cross_sections = parse_TOPBase_cross_sections(path, norad_format)

    nist_levels_path = joinpath(data_dir, "NIST_energy_levels", string(spec) * ".csv")
    nist_levels = if isfile(nist_levels_path)
        parse_NIST_energy_levels(nist_levels_path)
    else
        throw(ArgumentError(nist_levels_path * " doesn't exist"))
    end

    Z = get_atoms(spec)[1]
    ionization_energy = Korg.ionization_energies[Z][spec.charge+1]

    single_species_bf_cross_section(cross_sections, nist_levels, ionization_energy,
                                    Korg.default_partition_funcs[spec], λs, Ts)
end
function single_species_bf_cross_section(cross_sections, energy_levels, ionization_energy, U, λs,
                                         Ts)
    # convert λ_vals to photon energies in eV
    photon_energies = (hplanck_eV * c_cgs) ./ λs

    # precompute Temperature-dependent constant
    Us = U.(log.(Ts))
    β = 1.0 ./ (kboltz_eV .* Ts)

    #ionization_energy_topbase = maximum(E for (E, ) in values(cross_sections))
    #display("ionization energy: $(ionization_energy) (NIST) $(ionization_energy_topbase) (theoretical)")

    # prepare the output array where results will be accumulated.  This the cross section obtained 
    # by summing over electron states
    total_sigma = zeros(eltype(photon_energies), (length(photon_energies), length(Ts)))
    total_weight = zeros(length(Ts))

    for (state, (topbase_binding_energy, Es, σs)) in pairs(cross_sections)
        #the topbase binding energies are not precise.  Use NIST empirical energy levels instead.
        energy_level, g = if state in keys(energy_levels)
            energy_levels[state]
        else
            #println("skipping $state")
            continue
        end

        empirical_binding_energy = ionization_energy - energy_level #in eV

        # remove bb energies
        mask = Es .>= topbase_binding_energy
        Es = Es[mask]
        σs = σs[mask]

        # adjust Es (photon energies) to match empirical binding energy
        Es .+= empirical_binding_energy - topbase_binding_energy
        deduplicate_knots!(Es; move_knots=true) #shift repeat Es to the next float

        σ_itp = linear_interpolation(Es, σs; extrapolation_bc=0.0)

        # g*exp(-βε)/U at each temperature
        weights = g .* exp.(-energy_level .* β) ./ Us
        total_weight += weights

        # sum contributions from energy levels and adjust for stimulated recombination
        @. total_sigma += σ_itp(photon_energies) * (1.0 - exp(-photon_energies * β')) * weights'
    end

    #display(total_weight)

    total_sigma
end
