# note: this script should be run within a julia environment with these packes installed,
# not the Korg environment
using Korg, Downloads, CSV, DataFrames, HDF5

include("download_energy_levels.jl")

println("downloading energy levels from NIST...")
@time dfs = download_levels_from_NIST()

U(logT, df) = sum(@. df.g * exp(-df.level/(Korg.kboltz_eV * 10^logT)))
logTs = 0:0.01:5
h5write("partition_funcs.h5", "logT_min", logTs[1])
h5write("partition_funcs.h5", "logT_step", step(logTs))
h5write("partition_funcs.h5", "logT_max", logTs[end])

for (spec, df) in dfs
    Us = U.(logTs, Ref(df))
    h5write("partition_funcs.h5", string(spec), Us)
end

println("tabulated partition funcs written to partition_funcs.h5")