using Korg, Test, Logging, HDF5, ForwardDiff, FiniteDiff, TimerOutputs

@testset "Korg tests" begin

    # tools for testing: assert_allclose and assert_allclose_grid
    include("utilities.jl")

    # We use the TimerOutputs package to print the time and allocations for each top-level testset.
    # Inner testsets of particular interest can also be timed by wrapping them in a @timeit macro.

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
    @timeit "fit" include("fit.jl")                                   # slow
    @timeit "autodiff" include("autodiff.jl")                         # slow
    @timeit "autodiffable_conv" include("autodiffable_conv.jl")
    @timeit "atmosphere" include("atmosphere.jl")                     # slow
    @timeit "abundances" include("abundances.jl")
    @timeit "synthesize" include("synthesize.jl")
    @timeit "prune_linelist" include("prune_linelist.jl")
    @timeit "utils" include("utils.jl")
    @timeit "line_absorption" include("line_absorption.jl")
    @timeit "qfactors" include("qfactors.jl")
    print_timer()
end #top-level testset
