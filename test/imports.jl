using Korg, Test, Logging, HDF5, ForwardDiff, FiniteDiff, Aqua

println("running test suite with ", Threads.nthreads(), " threads.")

# tools for testing: assert_allclose and assert_allclose_grid
include("utilities.jl")
