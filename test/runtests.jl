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

@testset "linelist" begin
    @test SSSynth.parse_species_code("01.00") == "H_I"
    @test SSSynth.parse_species_code("02.01") == "He_II"

    @testset "linelist stub" begin
        linelist = SSSynth.read_line_list("data/gfallvac08oct17.stub.dat")
        @test length(linelist) == 999
        @test linelist[1].wl == 7232.0699
        @test linelist[1].log_gf == -0.826
        @test linelist[1].species == "Be_II"
        @test linelist[1].E == 140020.580
        @test linelist[1].log_gamma_rad == 7.93
        @test linelist[1].log_gamma_stark == -2.41
        @test linelist[1].log_gamma_vdW == -6.91
    end
end
