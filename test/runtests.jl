using Korg, Test, HDF5, ForwardDiff, FiniteDiff

@testset "Korg tests" begin

# tools for testing: assert_allclose and assert_allclose_grid
include("utilities.jl") 

# tests for specific parts of the code broken out into their own files. As you add tests, do it 
# this way.
include("molecular_cross_sections.jl")
include("cubic_splines.jl")
include("transfer.jl")
include("species.jl")
include("interval.jl")
include("continuum_absorption.jl") # test this after the "Interval" testset
include("partition_funcs.jl")
include("statmech.jl")
include("linelist.jl")
include("fit.jl")
include("autodiff.jl")
include("autodiffable_conv.jl")
include("atmosphere.jl")
include("synthesize.jl")
include("prune_linelist.jl")
include("utils.jl")
include("line_absorption.jl")


end #top-level testset
