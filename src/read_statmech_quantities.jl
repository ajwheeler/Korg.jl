using ..CubicSplines: CubicSpline

"""
    setup_ionization_energies([filename])

Parses the table of ionization energies and returns it as a dictionary mapping elements to
their ionization energies, `[χ₁, χ₂, χ₃]` in eV.
"""
function setup_ionization_energies(fname=joinpath(_data_dir, 
                                                  "BarklemCollet2016-ionization_energies.dat"))
    open(fname, "r") do f
        d = Dict{UInt8, Vector{Float64}}()
        for line in eachline(f)
            if line[1] != '#'        
                toks = split(strip(line))
                Z = parse(UInt8, toks[1])
                #the second token is the atomic symbol, which we ignore
                d[Z] = parse.(Float64, toks[3:end])
            end
        end
        d
    end
end

"""
    setup_partition_funcs()

Returns a Dict holding the default (Barklem & Collet 2016) partition functions.

The partition function is linearly interpolated from values tabulated for the paper. The partition 
functions can only be evaluated at temperatures between 1e-5 K and 1e4 K. The tabulated values were 
NOT computed with a cutoff principle quantum number. This can lead to a couple percent error for 
some neutral species in solar type stars (I think the problem becomes worse in hotter stars). See 
the paper for more details.

Note that the paper seems to suggest that we could actually use cubic interpolation. That should be
revisited in the future.
"""
function setup_partition_funcs()
    mol_path = joinpath(_data_dir, "BarklemCollet2016-molecular_partition.dat")
    merge(read_Barklem_Collet_table(mol_path), load_atomic_partition_functions())
end

"""
    setup_equilibrium_constants()

Returns a Dict holding the default (Barklem & Collet 2016) log equilibrium constants, which are in
pressure units.

In equation 7 in this paper, m is the reduced mass, m₁m₂/(m₁ + m₂).  This can be verified by 
re-deriving the equilibrium constants from the partition functions and dissolution constants 
provided by the paper.
"""
setup_equilibrium_constants() = read_Barklem_Collet_table(
                                 joinpath(_data_dir, "BarklemCollet2016-equilibrium_constants.dat"),
                                 transform=x->x+1) #convert from log(mks) to log(cgs)

"""
    functon read_Barklem_Collet_table([transform=identity])

Constructs a Dict holding tables containing partition function or equilibrium constant values across
ln(temperature).

The optional argument, `transform` is applied to each value. It is used to perform unit conversions,
when loading equilibrium constants
"""
function read_Barklem_Collet_table(fname; transform=identity)
    temperatures = Vector{Float64}()
    data_pairs = Vector{Tuple{Species,Vector{Float64}}}()
    open(fname, "r") do f
        for line in eachline(f)
            if (length(line)>=9) && contains(line, "T [K]")
                append!(temperatures, parse.(Float64,split(strip(line[10:length(line)]))))
            elseif line[1] == '#'
                continue
            else # add entries to the dictionary
                row_entries = split(strip(line))
                species_code = popfirst!(row_entries)
                #ignore deuterium and ionized molecules
                if species_code[1:2] != "D_" && species_code[end] != '+' && species_code[end] != '-' 
                    push!(data_pairs, (Species(species_code), 
                                       transform.(parse.(Float64, row_entries))))
                end
            end
        end
    end

    map(data_pairs) do (species, vals)
        species, CubicSpline(log.(temperatures), vals)
    end |> Dict
end


"""
    load_atomic_partition_functions()

Loads saved tabulated values for atomic partition functions from disk. Returns a dictionary mapping
species to interpolators over log(T).
"""
function load_atomic_partition_functions(filename=joinpath(_data_dir, "atomic_partition_funcs", 
                                         "partition_funcs.h5"))
    partition_funcs = Dict{Species, Any}()

    logT_min = h5read(filename, "logT_min")
    logT_step = h5read(filename, "logT_step")
    logT_max = h5read(filename, "logT_max")
    logTs = logT_min : logT_step : logT_max

    for elem in Korg.atomic_symbols, ionization in ["I", "II", "III"]
        if (elem == "H" && ionization != "I") || (elem == "He" && ionization == "III")
            continue
        end
        spec = elem*" "*ionization
        partition_funcs[Species(spec)] = CubicSpline(logTs, h5read(filename, spec))
    end

    #handle the cases with bare nuclei
    all_ones = ones(length(logTs))
    partition_funcs[species"H II"] = CubicSpline(logTs, all_ones)
    partition_funcs[species"He III"] = CubicSpline(logTs, all_ones)

    partition_funcs
end