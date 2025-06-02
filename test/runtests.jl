using Korg, Test, Logging, HDF5, ForwardDiff, FiniteDiff, Aqua

println("running test suite with ", Threads.nthreads(), " threads.")

@testset "Korg tests" verbose=true showtiming=true begin

    # tools for testing: assert_allclose and assert_allclose_grid
    include("utilities.jl")

    include("wavelengths.jl")
    include("molecular_cross_sections.jl")
    include("cubic_splines.jl")
    include("transfer.jl")
    include("species.jl")
    include("interval.jl")
    include("continuum_absorption.jl") # test this after the "Interval" testset
    include("statmech.jl")
    include("linelist.jl")
    include("fit.jl")                            # slow
    include("autodiff.jl")                       # slow
    include("autodiffable_conv.jl")
    include("atmosphere.jl")                     # slow
    include("abundances.jl")
    include("synthesize.jl")
    include("synth.jl")
    include("prune_linelist.jl")
    include("utils.jl")
    include("line_absorption.jl")
    include("qfactors.jl")
    @testset "Aqua automated checks" begin
        # see https://github.com/JuliaTesting/Aqua.jl/issues/77 for why I'm doing it this way,
        # (basically the default bahavior is different and this is what we want to avoid errors in
        # deps that we can't fix.)
        # unbound args has a false positive for Korg.Line because it's using a heuristic.  That
        # constructor should probably be less convoluted anyway, but for now skip those checks.
        Aqua.test_all(Korg; ambiguities=false, unbound_args=false)
        Aqua.test_ambiguities(Korg)
    end
end #top-level testset
