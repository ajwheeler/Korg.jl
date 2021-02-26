using Interpolations

_data_dir = joinpath(@__DIR__, "../data") # may want to pull this out of here later

# We had the option of downloading the file as a FITS table, but when we do that, we lose the
# information about the temperatures. It might make sense to repackage the data in the future.
function _load_atomic_data(fname)
    temperatures = Vector{Float64}()
    data_pairs = Vector{Tuple{AbstractString,Vector{Float64}}}()

    open(fname, "r") do f
        for line in eachline(f)
            if (length(line)>=9) && (line[1:9] == "#   T [K]")
                append!(temperatures, map(x -> parse(Float64, x),
                                          split(strip(line[10:length(line)]))))
            elseif line[1] == '#'
                continue
            else # add entries to the dictionary
                row_entries = split(strip(line))
                species = popfirst!(row_entries)
                push!(data_pairs,
                      (species, Vector(map(x -> parse(Float64, x), row_entries))))
            end
        end
    end
    return (temperatures, data_pairs)
end

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
function setup_atomic_partition_funcs()
    fname = joinpath(_data_dir, "BarklemCollet2016-atomic_partition.dat")
    temperatures, data_pairs = _load_atomic_data(fname)

    # note that elem is a struct that unpacks to (key, partition_function_values)
    _helper = elem -> (elem[1],
                       Interpolations.LinearInterpolation(temperatures, elem[2],
                                                          extrapolation_bc=Interpolations.Throw()))
    Dict(map(_helper,data_pairs))
end
