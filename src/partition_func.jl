using Interpolations: LinearInterpolation, Throw

"""
    setup_partition_funcs()

Returns a Dict holding the default (Barklem & Collet 2016) partition functions.
"""
function setup_partition_funcs(atoms=joinpath(_data_dir, "BarklemCollet2016-atomic_partition.dat"),
                              mols=joinpath(_data_dir, "BarklemCollet2016-molecular_partition.dat"))
    atomUs, molUs = read_partition_funcs.([atoms, mols])
    new_molUs = typeof(molUs)()
    #use only neutral molecules, add "I"
    for (m, f) in molUs
        if !('+' in m) && !('-' in m)
            new_molUs[m*"_I"] = f
        end
    end
    merge(atomUs, new_molUs)
end

"""
    setup_equilibrium_constants()

Returns a Dict holding the default (Barklem & Collet 2016) log equilibrium constants.

There are two things of note about about equation 7 in this paper.  The first is that it defines K 
as the reciprocal of what is actually supplied in the data table.  The second is that in order for 
it to be correct, m must be interpreted as a reduced mass, m₁m₂/(m₁ + m₂).  This can be verified by 
re-deriving the equilibrium constants from the partition functions and dissolution constants 
provided by the paper.
"""
setup_equilibrium_constants() = read_partition_funcs(
                                 joinpath(_data_dir, "BarklemCollet2016-equilibrium_constants.dat"))

"""
    functon setup_atomic_partition_funcs()

Constructs a Dict holding functions for various atomic species that computes the partition function
at an input temperature.

In more detail, the partition function is linearly interpolated from values tabulated from Barklem 
& Collet (2016). The partition functions can only be evaluated at temperatures between 1e-5 K and 
1e4 K. The tabulated values were NOT computed with a cutoff principle quantum number. This can lead
to a couple percent error for some neutral species in solar type stars (I think the problem becomes
worse in hotter stars). See the paper for more details.

Note that the paper seems to suggest that we could actually use cubic interpolation. That should be
revisited in the future.
"""
function read_partition_funcs(fname)
    temperatures = Vector{Float64}()
    data_pairs = Vector{Tuple{AbstractString,Vector{Float64}}}()
    open(fname, "r") do f
        for line in eachline(f)
            if (length(line)>=9) && contains(line, "T [K]")
                append!(temperatures, parse.(Float64,split(strip(line[10:length(line)]))))
            elseif line[1] == '#'
                continue
            else # add entries to the dictionary
                row_entries = split(strip(line))
                species = popfirst!(row_entries)
                if species[end] != '+' #ignore ionized molecules
                    push!(data_pairs, (species, parse.(Float64, row_entries)))
                end
            end
        end
    end

    map(data_pairs) do (species, vals)
        species, LinearInterpolation(temperatures, vals, extrapolation_bc=Throw())
    end |> Dict
end
