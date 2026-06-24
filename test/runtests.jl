# This file is automatically run by Pkg.test() when the Korg environment is activated.
# This is the standard way that tests work in Julia packages. It's how the CI works, and it's a way
# to run the whole test suite locally.
#
# When testing locally, another recommended workflow is open a REPL and run
#     cd("test/"); import Pkg; Pkg.activate("."); using Revise; include("imports.jl")
# Then include individual test files as needed to run them. You can iterate on
# Korg and on the tests, rerunning as needed with minimal recompilation

# Imports are in a separate file to facilitate running individual test files.
include("imports.jl")

@testset "Korg tests" verbose=true showtiming=true begin
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
    include("newton.jl")
    include("atmosphere.jl")                     # slow
    include("abundances.jl")
    include("synthesize.jl")
    include("synth.jl")
    include("prune_linelist.jl")
    include("utils.jl")
    include("line_absorption.jl")
    include("qfactors.jl")
    @testset "Aqua automated checks" begin
        # unbound args has a false positive for Korg.Line because it's using a heuristic.  That
        # constructor should probably be less convoluted anyway, but for now skip those checks.
        Aqua.test_all(Korg; unbound_args=false)
    end
end #top-level testset
