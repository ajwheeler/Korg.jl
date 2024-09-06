# note: this script should be run within a julia environment with these packes installed,
# not the Korg environment
using Korg, Downloads, CSV, DataFrames, HDF5

include("download_energy_levels.jl")

println("downloading energy levels from NIST...")
@time dfs = download_levels_from_NIST()

U(lnT, df) = sum(@. df.g * exp(-df.level / (Korg.kboltz_eV * exp(lnT))))

lnTs = (0:0.025:5) * log(10) # natural log of 1 to 10,000

h5write("partition_funcs.h5", "logT_min", lnTs[1])
h5write("partition_funcs.h5", "logT_step", step(lnTs))
h5write("partition_funcs.h5", "logT_max", lnTs[end])

for (spec, df) in dfs
    Us = U.(lnTs, Ref(df))
    h5write("partition_funcs.h5", string(spec), Us)
end

println("tabulated partition funcs written to partition_funcs.h5")
