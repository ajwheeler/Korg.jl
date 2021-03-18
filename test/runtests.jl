using SSSynth
using Test

@testset "ionization energies" begin
    @test length(SSSynth.ionization_energies) == 92
    @test SSSynth.ionization_energies["H"] == [13.5984, -1.000, -1.000]
    @test SSSynth.ionization_energies["Ru"] == [7.3605, 16.760, 28.470]
    @test SSSynth.ionization_energies["U"] == [6.1940, 11.590, 19.800]
end 

@testset "saha" begin s = [SSSynth.saha(SSSynth.ionization_energies["N"], 
                      [SSSynth.partition_funcs["N_I"], 
                       SSSynth.partition_funcs["N_II"], 
                       SSSynth.partition_funcs["N_III"]], 
                      T, 1.0) 
         for T in 1:100:10000]
    @test issorted(first.(s), rev=true)
    @test issorted(last.(s))
    @test issorted(s[1], rev=true)
end

@testset "lines" begin
    @testset "species codes" begin
        @test SSSynth.parse_species_code("01.00") == "H_I"
        @test SSSynth.parse_species_code("02.01") == "He_II"
    end

    linelist = SSSynth.read_line_list("data/gfallvac08oct17.stub.dat")
    @testset "linelist parsing" begin
        @test length(linelist) == 988
        @test linelist[1].wl ≈ 72320.699
        @test linelist[1].log_gf == -0.826
        @test linelist[1].species == "Be_II"
        @test linelist[1].wavenumber == 140020.580
        @test linelist[1].log_gamma_rad == 7.93
        @test linelist[1].log_gamma_stark == -2.41
        @test linelist[1].log_gamma_vdW == -6.91
    end

    @testset "line profile" begin
        @test issorted(SSSynth.line_profile(5000.0, SSSynth.atomic_masses["Be"], linelist[1], 
                                            7230 : 0.01 : linelist[1].wl))
        @test issorted(SSSynth.line_profile(5000.0, SSSynth.atomic_masses["Be"], linelist[1], 
                                            linelist[1].wl : 0.01 : 7235), rev=true)

        s =  sum(SSSynth.line_profile(5000.0, SSSynth.atomic_masses["Be"], linelist[1],
                                      72300 : 0.1 : 72350))/10
        @test 0.999 < s <= 1.0
    end
end

@testset "atmosphere" begin
    #the MARCS solar model atmosphere
    atmosphere = SSSynth.read_model_atmosphere("data/sun.krz")
    @test length(atmosphere) == 56
    @test issorted(first.(atmosphere))
    @test atmosphere[1].colmass == 9.747804143e-3
    @test atmosphere[1].temp == 4066.8
    @test atmosphere[1].electron_density == 3.76980e10
    @test atmosphere[1].number_density == 4.75478e14
    @test atmosphere[1].density == 1.00062e-9
end

include("continuum_opacity.jl")
