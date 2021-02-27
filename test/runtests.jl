using SSSynth
using Test

@testset "ionization energies" begin
    ionEs = SSSynth.setup_ionization_energies()
    @test length(ionEs) == 92
    @test ionEs["H"] == [13.5984, -1.000, -1.000]
    @test ionEs["Ru"] == [7.3605, 16.760, 28.470]
    @test ionEs["U"] == [6.1940, 11.590, 19.800]
end 
