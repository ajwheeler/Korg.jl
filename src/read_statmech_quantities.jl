using ..CubicSplines: CubicSpline
using HDF5, CSV

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
    setup_partition_funcs_and_equilibrium_constants()

Returns two dictionaries. One holding the default patition functions, and one holding the default
log10 equilibrium constants.

# Default partition functions
The partition functions are custom (calculated from NIST levels) for atoms, from Barklem & 
Collet 2016 for diatomic molecules, and from [exomol](https://exomol.com) for polyatomic molecules.
For each molecule, we include only the most abundant isotopologue.

Note than none of these parition functions include plasma effects, e.g. via the Mihalas Hummer 
Daeppen occupation probability formalism. They are for isolated species.
This can lead to a couple percent error for neutral alkalis and to greater errors for hydrogen in 
some atmospheres, particularly those of hot stars.

# Default equilibrium constants
Molecules have equilibrium constants in addition to patition functions.  For the diatomics, these 
are provided by Barklem and Collet, which extensively discusses the dissociation energies.  For 
polyatomics, we calculate these ourselves, using atomization energies calculated from the enthalpies 
of formation at 0K from [NIST's CCCDB](https://cccbdb.nist.gov/hf0k.asp).

Korg's equilbrum constants are in terms of partial pressures, since that's what Barklem and Collet 
provide.
"""
function setup_partition_funcs_and_equilibrium_constants()
    mol_U_path = joinpath(_data_dir, "BarklemCollet2016-molecular_partition.dat")
    partition_funcs = merge(read_Barklem_Collet_table(mol_U_path), 
                           load_atomic_partition_functions(),
                           load_exomol_partition_functions())

    BC_Ks = read_Barklem_Collet_table(joinpath(_data_dir, "BarklemCollet2016-equilibrium_constants.dat"),
                                      transform=x->x+1) #convert from log(mks) to log(cgs)

    # equlibrium constants for polyatomics (done here because the partition functions are used)
    atomization_Es = CSV.File(joinpath(_data_dir, "polyatomic_partition_funcs", "atomization_energies.csv"))
    polyatomic_Ks = map(zip(Korg.Species.(atomization_Es.spec), atomization_Es.energy)) do (spec, D00)
        D00 *= 0.01036 # convert from kJ/mol to eV
        # this let block slightly improves performance. 
        # https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured
        calculate_logK = let parition_funcs=partition_funcs
            function logK(logT)
                Zs = get_atoms(spec)
                log_Us_ratio = log10(prod([partition_funcs[Species(Formula(Z), 0)](logT) for Z in Zs])
                                     / partition_funcs[spec](logT))
                log_masses_ratio = sum([log10(atomic_masses[Z]) for Z in Zs]) - log10(get_mass(spec))
                T = exp(logT)
                log_translational_U_factor = 1.5*log10(2π*kboltz_cgs*T/hplanck_cgs^2)
                # this is log number-density equilbrium constant
                log_nK = ((length(Zs)-1)*log_translational_U_factor 
                          + 1.5*log_masses_ratio + log_Us_ratio - D00/(kboltz_eV*T*log(10)))
                # compute the log of the partial-pressure equilibrium constant, log10(pK)
                log_nK + (length(Zs)-1)*log10(kboltz_cgs*T)
            end
        end
        spec => calculate_logK
    end |> Dict

    # We do this rather than merge(BC_Ks, polyatomic_Ks) because the latter would produce a valtype
    # of Any.
    equilibrium_constants = Dict{Korg.Species, Union{valtype(BC_Ks), valtype(polyatomic_Ks)}}()
    merge!(equilibrium_constants, BC_Ks)
    merge!(equilibrium_constants, polyatomic_Ks)
    

    partition_funcs, equilibrium_constants
end

"""
    functon read_Barklem_Collet_table(;transform=identity)

Constructs a Dict holding tables containing partition function or equilibrium constant values across
ln(temperature).  Applies transform (which you can use to, e.g. change units) to each example.
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
                # ignore deuterium and ionized molecules
                if species_code[1:2] != "D_" && species_code[end] != '+' && species_code[end] != '-' 
                    push!(data_pairs, (Species(species_code), 
                                       transform.(parse.(Float64, row_entries))))
                end
            end
        end
    end

    map(data_pairs) do (species, vals)
        # We want to extrapolate the molecular partition funcs because they are only defined up to 
        # 10,000 K. There won't be any molecules past that temperature anyway.
        species, CubicSpline(log.(temperatures), vals; extrapolate=true)
    end |> Dict
end


"""
    load_atomic_partition_functions()

Loads saved tabulated values for atomic partition functions from disk. Returns a dictionary mapping
species to interpolators over log(T).
"""
function load_atomic_partition_functions(filename=joinpath(_data_dir, "atomic_partition_funcs", 
                                         "partition_funcs.h5"))
    partition_funcs = Dict{Species, CubicSplines.CubicSpline{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64}}()

    logT_min = h5read(filename, "logT_min")
    logT_step = h5read(filename, "logT_step")
    logT_max = h5read(filename, "logT_max")
    logTs = collect(logT_min : logT_step : logT_max)

    for elem in Korg.atomic_symbols, ionization in ["I", "II", "III"]
        if (elem == "H" && ionization != "I") || (elem == "He" && ionization == "III")
            continue
        end
        spec = elem*" "*ionization
        # this is flat "extrapolation", not linear or cubic
        partition_funcs[Species(spec)] = CubicSpline(logTs, h5read(filename, spec))
    end

    #handle the cases with bare nuclei
    all_ones = ones(length(logTs))
    partition_funcs[species"H II"] = CubicSpline(logTs, all_ones)
    partition_funcs[species"He III"] = CubicSpline(logTs, all_ones)

    partition_funcs
end

"""
    load_exomol_partition_functions()

Loads the exomol partition funcstions for polyatomic molecules from the HDF5 archive. Returns a 
dictionary mapping species to interpolators over log(T).
"""
function load_exomol_partition_functions()
    h5open(joinpath(_data_dir, "polyatomic_partition_funcs", "polyatomic_partition_funcs.h5")) do f
        map(f) do group
            spec = Korg.Species(HDF5.name(group)[2:end]) # cut off leading '/'
            Ts, Us = read(group["temp"]), read(group["partition_function"])
            spec, CubicSpline(log.(Ts), Us, extrapolate=true)
        end |> Dict
    end
end

#load data when the package is imported. 
const ionization_energies = setup_ionization_energies()
const default_partition_funcs, default_log_equilibrium_constants = setup_partition_funcs_and_equilibrium_constants()