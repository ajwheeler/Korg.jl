using BenchmarkTools, Korg

const SUITE = BenchmarkGroup()

include("synthesis.jl")
include("fit.jl")
include("chemical_equilibrium.jl")

results = run(SUITE)
println(results)
