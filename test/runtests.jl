using SSSynth
using Test

@testset "ionization energies" begin
    @test length(SSSynth.ionization_energies) == 92
    @test SSSynth.ionization_energies["H"] == [13.5984, -1.000, -1.000]
    @test SSSynth.ionization_energies["Ru"] == [7.3605, 16.760, 28.470]
    @test SSSynth.ionization_energies["U"] == [6.1940, 11.590, 19.800]
end


"""
Compute nₑ (number density of free electrons) in a pure Hydrogen atmosphere, where `nH_tot` is the
total number density of H I and H II (in cm⁻³), the temperature is `T`, and `HI_partition_val` is
the the value of the H I partition function.

This is a relatively naive implementation. More numerically stable solutions exist.
"""
function electron_ndens_Hplasma(nH_tot, T, H_I_partition_val = 2.0)
    # Define the Saha equation as: nₑ*n_{H II} / n_{H I} = RHS
    # coef ∼ 4.829e15
    coef = 2.0 * (2.0*π*SSSynth.electron_mass_cgs*SSSynth.kboltz_cgs / SSSynth.hplanck_cgs^2)^1.5
    RHS = coef * T^1.5 * exp(-SSSynth.RydbergH_eV/(SSSynth.kboltz_eV*T))/H_I_partition_val
    # In a pure Hydrogen atmosphere: nₑ = n_{H II}. The Saha eqn becomes:  nₑ²/(nH_tot - ne) = RHS
    # We recast the Saha eqn as: a*nₑ² + b*nₑ + c = 0 and compute the coefficients
    a, b, c = (1.0, RHS, -1*RHS*nH_tot)
    # solve quadratic equation. Since b is always positive and c is always negative:
    #    (-b + sqrt(b²-4*a*c))/(2*a) is always ≥ 0
    #    (-b - sqrt(b²-4*a*c))/(2*a) is always negative
    nₑ = (-b + sqrt(b*b-4*a*c))/(2*a)
    nₑ
end

@testset "saha" begin
    @testset "pure Hydrogen atmosphere" begin
        nH_tot = 1e15
        # specify χs and Us to decouple this testset from other parts of the code
        χs, Us = [SSSynth.RydbergH_eV, -1.0], [T -> 2.0, T -> 1.0]
        # iterate from less than 1% ionized to more than 99% ionized
        for T in [3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 1e4, 1.1e4, 1.2e4, 1.3e4, 1.4e4, 1.5e5]
            nₑ = electron_ndens_Hplasma(nH_tot, T, 2.0)
            weights = SSSynth.saha(χs, Us, T, nₑ)
            @test length(weights) == 2
            @test (weights[1] + weights[2]) ≈ 1.0 rtol = 1e-15
            rtol = (T == 1.5e5) ? 1e-9 : 1e-14
            @test weights[2] ≈ (nₑ/nH_tot) rtol = rtol
        end
    end

    @testset "monotonic N_I and N_III Temperature dependence" begin
        s = [SSSynth.saha(SSSynth.ionization_energies["N"], [SSSynth.partition_funcs["N_I"],
                                                             SSSynth.partition_funcs["N_II"],
                                                             SSSynth.partition_funcs["N_III"]],
                          T, 1.0) for T in 1:100:10000]
        @test issorted(first.(s), rev=true)
        @test issorted(last.(s))
        @test issorted(s[1], rev=true)
    end
end

@testset "lines" begin
    @testset "species codes" begin
        @test SSSynth.parse_species_code("01.00") == "H_I"
        @test SSSynth.parse_species_code("02.01") == "He_II"
    end

    linelist = SSSynth.read_line_list("data/gfallvac08oct17.stub.dat")
    @testset "linelist parsing" begin
        @test length(linelist) == 999
        @test linelist[1].wl ≈ 72320.699
        @test linelist[1].log_gf == -0.826
        @test linelist[1].species == "Be_II"
        @test linelist[1].E == 140020.580
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
    @test atmosphere[1].tau == 9.747804143e-3
    @test atmosphere[1].temp == 4066.8
    @test atmosphere[1].electron_density == 3.76980e10
    @test atmosphere[1].number_density == 4.75478e14
    @test atmosphere[1].density == 1.00062e-9
end

include("continuum_opacity.jl")
