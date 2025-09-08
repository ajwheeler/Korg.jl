using BenchmarkTools, Korg

const SUITE = BenchmarkGroup()

include("synthesis.jl")
include("fit.jl")
include("chemical_equilibrium.jl")

@show SUITE

results = run(SUITE)
println(results)
