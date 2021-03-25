using SSSynth
using Test

include("continuum_opacity.jl")

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
    @testset "line lists" begin 
        @testset "species codes" begin
            @test SSSynth.parse_species_code("01.00") == "H_I"
            @test SSSynth.parse_species_code("02.01") == "He_II"
        end

        @test_throws ArgumentError SSSynth.read_line_list(""; format="abc")

        kurucz_linelist = SSSynth.read_line_list("data/gfallvac08oct17.stub.dat")
        @testset "kurucz linelist parsing" begin
            @test all((l->l.E_upper > l.E_lower).(kurucz_linelist))
            @test issorted(kurucz_linelist, by=l->l.wl)
            @test length(kurucz_linelist) == 988
            @test kurucz_linelist[1].wl ≈ 72320.699
            @test kurucz_linelist[1].log_gf == -0.826
            @test kurucz_linelist[1].species == "Be_II"
            @test kurucz_linelist[1].E_upper ≈ 17.531776
            @test kurucz_linelist[1].E_lower ≈ 17.360339371573698
            @test kurucz_linelist[1].log_gamma_rad == 7.93
            @test kurucz_linelist[1].log_gamma_stark == -2.41
            @test kurucz_linelist[1].log_gamma_vdW == -6.91
        end

        vald_linelist = SSSynth.read_line_list("data/Ylines.vald"; format="vald")
        @testset "vald linelist parsing" begin
            @test all((l->l.E_upper > l.E_lower).(vald_linelist))
            @test issorted(vald_linelist, by=l->l.wl)
            @test length(vald_linelist) == 4584
            @test vald_linelist[1].wl ≈ 3002.20106
            @test vald_linelist[1].log_gf == -1.132
            @test vald_linelist[1].species == "Y_II"
            @test vald_linelist[1].E_lower ≈ 3.3757
            @test vald_linelist[1].E_upper ≈ 7.5055 
            @test vald_linelist[1].log_gamma_rad == 8.620
            @test vald_linelist[1].log_gamma_stark == -5.580
            @test vald_linelist[1].log_gamma_vdW == -7.710
        end

        @test typeof(vald_linelist) == typeof(kurucz_linelist)
    end
        


    @testset "line profile" begin
        linelist = SSSynth.read_line_list("data/gfallvac08oct17.stub.dat")

        @test issorted(SSSynth.line_profile(5000.0, SSSynth.atomic_masses["Be"], linelist[1], 
                                            7230 : 0.01 : linelist[1].wl))
        @test issorted(SSSynth.line_profile(5000.0, SSSynth.atomic_masses["Be"], linelist[1], 
                                            linelist[1].wl : 0.01 : 7235), rev=true)

        s =  sum(SSSynth.line_profile(5000.0, SSSynth.atomic_masses["Be"], linelist[1],
                                      72300 : 0.1 : 72350))/10
        #must convert from cm^-1 to Å^-1
        @test 0.999 < s * 1e-8 <= 1.0
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

@testset "synthesis" begin

    @testset "calculate absolute abundances" begin
        @test_throws ArgumentError SSSynth.get_absolute_abundances(["H"], 0.0, Dict("H"=>13))

        @testset for metallicity in [0.0, 1.0], A_X in [Dict(), Dict("C"=>9)]
            for elements in [SSSynth.atomic_symbols, ["H", "He", "C", "Ba"]]
                nxnt = SSSynth.get_absolute_abundances(elements, metallicity, A_X)

                #abundances for the right set of elementns
                @test Set(elements) == Set(keys(nxnt))

                #correct absolute abundances?
                if "C" in keys(A_X)
                    @test log10(nxnt["C"]/nxnt["H"]) + 12 ≈ 9
                end
                @test log10(nxnt["He"]/nxnt["H"]) + 12 ≈ SSSynth.solar_abundances["He"]
                @test log10(nxnt["Ba"]/nxnt["H"]) + 12 ≈ 
                    SSSynth.solar_abundances["Ba"] + metallicity

                #normalized?
                if elements == SSSynth.atomic_symbols
                    @test sum(values(nxnt)) ≈ 1
                else
                    @test sum(values(nxnt)) < 1
                end
            end
        end
    end

    @testset "number densities" begin
        abundances = Dict("C"=>1e-8, "H"=>0.9)
        nₜ = 1e15
        n = SSSynth.per_species_number_density(nₜ, nₜ*1e-3, 4500.0, abundances)
        @test Set(keys(n)) == Set(["C_I", "C_II", "C_III", "H_I", "H_II"])
        @test n["C_III"] < n["C_II"] < n["C_I"] < n["H_II"] < n["H_I"]
        @test n["C_III"] + n["C_II"] + n["C_I"] ≈ abundances["C"] * nₜ
        @test n["H_II"] + n["H_I"] ≈ abundances["H"] * nₜ
    end

    @testset "trapezoid rule" begin
        pdf(x) = exp(-1/2 * x^2) / sqrt(2π)
        xs = -10:0.1:10
        @test SSSynth.trapezoid_rule(xs, pdf.(xs) * 0.1) - 1.0 < 1e-5
    end
    
end
