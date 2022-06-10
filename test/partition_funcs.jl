@testset "O I-III and CN partition functions are (nearly) monotonic in T" begin
    lnTs = 0:0.1:log(100_000.0)
    nearly_monotonic(Us) = all(diff(Us) .> -1e-4)
    @test nearly_monotonic(Korg.partition_funcs[Korg.species"O_I"].(lnTs))
    @test nearly_monotonic(Korg.partition_funcs[Korg.species"O_II"].(lnTs))
    @test nearly_monotonic(Korg.partition_funcs[Korg.species"O_III"].(lnTs))
    @test nearly_monotonic(Korg.partition_funcs[Korg.species"CN_I"].(lnTs))
end

