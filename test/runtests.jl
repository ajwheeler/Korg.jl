using Korg, Test, HDF5, ForwardDiff, FiniteDiff, TimerOutputs

@testset "Korg tests" begin

# tools for testing: assert_allclose and assert_allclose_grid
include("utilities.jl") 

reset_timer!()
@timeit "molecular_cross_sections" include("molecular_cross_sections.jl")
@timeit "cubic_splines" include("cubic_splines.jl")
@timeit "transfer" include("transfer.jl")
@timeit "species" include("species.jl")
@timeit "interval" include("interval.jl")
@timeit "continuum_absorption" include("continuum_absorption.jl") # test this after the "Interval" testset
@timeit "partition_funcs" include("partition_funcs.jl")
@timeit "statmech" include("statmech.jl")
@timeit "linelist" include("linelist.jl")
@timeit "fit" include("fit.jl") # slow
@timeit "autodiff" include("autodiff.jl") # slow
@timeit "autodiffable_conv" include("autodiffable_conv.jl")
@timeit "atmosphere" include("atmosphere.jl")
@timeit "synthesize" include("synthesize.jl")
@timeit "prune_linelist" include("prune_linelist.jl")
@timeit "utils" include("utils.jl")
@timeit "line_absorption" include("line_absorption.jl")
print_timer()


end #top-level testset
