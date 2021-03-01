using SSSynth
using Test

@testset "ionization energies" begin
    @test length(SSSynth.ionization_energies) == 92
    @test SSSynth.ionization_energies["H"] == [13.5984, -1.000, -1.000]
    @test SSSynth.ionization_energies["Ru"] == [7.3605, 16.760, 28.470]
    @test SSSynth.ionization_energies["U"] == [6.1940, 11.590, 19.800]
end 

@testset "saha" begin
    s = [SSSynth.saha(SSSynth.ionization_energies["N"], 
                      [SSSynth.partition_funcs["N_I"], 
                       SSSynth.partition_funcs["N_II"], 
                       SSSynth.partition_funcs["N_III"]], 
                      T, 1.0) 
         for T in 1:100:10000]
    @test issorted(first.(s), rev=true)
    @test issorted(last.(s))
    @test issorted(s[1], rev=true)
end
