using SSSynth
using Test

include("continuum_opacity.jl")

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

@testset "O I-III and CN partition functions are monotonic in T" begin
    Ts = 1:100:10000
    @test issorted(SSSynth.partition_funcs["O_I"].(Ts))
    @test issorted(SSSynth.partition_funcs["O_II"].(Ts))
    @test issorted(SSSynth.partition_funcs["O_III"].(Ts))
    @test issorted(SSSynth.partition_funcs["CN_I"].(Ts))
end

@testset "stat mech" begin
    @testset "pure Hydrogen atmosphere" begin
        nH_tot = 1e15
        # specify χs and Us to decouple this testset from other parts of the code
        χs = Dict("H"=>[SSSynth.RydbergH_eV, -1.0, -1.0])
        Us = Dict(["H_I"=>(T -> 2.0), "H_II"=>(T -> 1.0)])
        # iterate from less than 1% ionized to more than 99% ionized
        for T in [3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 1e4, 1.1e4, 1.2e4, 1.3e4, 1.4e4, 1.5e5]
            nₑ = electron_ndens_Hplasma(nH_tot, T, 2.0)
            wII, wIII = SSSynth.saha_ion_weights(T, nₑ, "H", χs, Us)
            @test wIII ≈ 0.0 rtol = 1e-15
            rtol = (T == 1.5e5) ? 1e-9 : 1e-14
            @test wII/(1 + wII + wIII) ≈ (nₑ/nH_tot) rtol= rtol
        end
    end

    @testset "monotonic N ions Temperature dependence" begin
        weights = [SSSynth.saha_ion_weights(T, 1.0, "N", SSSynth.ionization_energies, 
                                            SSSynth.partition_funcs) for T in 1:100:10000]
        #N II + NIII grows with T === N I shrinks with T
        @test issorted(first.(weights) + last.(weights))
        
        # NIII grows with T
        @test issorted(last.(weights))
    end

    @testset "molecular equilibrium" begin
        #solar abundances
        abundances = SSSynth.get_absolute_abundances(SSSynth.atomic_symbols, 0.0, Dict())
        nₜ = 1e15 
        nₑ = 1e-3 * nₜ #arbitrary

        MEQs = SSSynth.molecular_equilibrium_equations(abundances, SSSynth.ionization_energies, 
                                                       SSSynth.partition_funcs, 
                                                       SSSynth.equilibrium_constants)

        #this should hold for the default atomic/molecular data
        @test Set(MEQs.atoms) == Set(SSSynth.atomic_symbols)

        n = SSSynth.molecular_equilibrium(MEQs, 5700.0, nₜ, nₑ)
        #make sure number densities are sensible
        @test n["C_III"] < n["C_II"] < n["C_I"] < n["H_II"] < n["H_I"]

        #total number of carbons is correct
        total_C = map(collect(keys(n))) do species
            println(species)
            if SSSynth.strip_ionization(species) == "C2"
                n[species] * 2
            elseif ((SSSynth.strip_ionization(species) == "C") || 
                    (SSSynth.ismolecule(species) && ("C" in SSSynth.get_atoms(species))))
                n[species]
            else
                0.0
            end
        end |> sum
        @test total_C ≈ abundances["C"] * nₜ
    end
end

@testset "lines" begin
    @testset "line lists" begin 
        @testset "species codes" begin
            @test SSSynth.parse_species_code("01.00") == "H_I"
            @test SSSynth.parse_species_code("02.01") == "He_II"
            @test SSSynth.parse_species_code("0608") == "CO"
        end

        @testset "strip ionization info" begin
            @test SSSynth.strip_ionization("H_I") == "H"
            @test SSSynth.strip_ionization("H_II") == "H"
            @test SSSynth.strip_ionization("CO") == "CO"
        end

        @testset "distinguish atoms from molecules" begin
            @test !SSSynth.ismolecule(SSSynth.strip_ionization("H_I"))
            @test !SSSynth.ismolecule(SSSynth.strip_ionization("H_II"))
            @test SSSynth.ismolecule(SSSynth.strip_ionization("CO"))

            @test !SSSynth.ismolecule("H_I")
            @test !SSSynth.ismolecule("H_II")
            @test SSSynth.ismolecule("CO")
        end

        @testset "break molecules into atoms" begin
            @test SSSynth.get_atoms("CO") == ("C", "O")
            @test SSSynth.get_atoms("C2") == ("C", "C")
            @test SSSynth.get_atoms("MgO") == ("Mg", "O")
            #nonsensical but it doesn't matter
            @test SSSynth.get_atoms("OMg") == ("O", "Mg")
            @test_throws ArgumentError SSSynth.get_atoms("hello world")
        end

        @test_throws ArgumentError SSSynth.read_line_list(""; format="abc")

        kurucz_linelist = SSSynth.read_line_list("data/gfallvac08oct17.stub.dat")
        @testset "kurucz linelist parsing" begin
            @test all((l->l.E_upper > l.E_lower).(kurucz_linelist))
            @test issorted(kurucz_linelist, by=l->l.wl)
            @test length(kurucz_linelist) == 988
            @test kurucz_linelist[1].wl ≈ 72320.699 * 1e-8
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
            @test vald_linelist[1].wl ≈ 3002.20106 * 1e-8
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

        @test issorted(SSSynth.line_profile(5000.0, SSSynth.atomic_masses["Be"], 1e5, linelist[1], 
                                            72.30e-5 : 1e-10 : linelist[1].wl))
        @test issorted(SSSynth.line_profile(5000.0, SSSynth.atomic_masses["Be"], 0.0, linelist[1], 
                                            linelist[1].wl : 1e-10 : 72.35e-5), rev=true)

        Δ = 1e-9 #cm (== 0.1 Å)
        s =  sum(SSSynth.line_profile(5000.0, SSSynth.atomic_masses["Be"], 2e5, linelist[1],
                                      72.300e-5 : Δ : 72.350e-5))
        #must convert from cm^-1 to Å^-1
        @test 0.999 < s * Δ <= 1.0
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

    @testset "trapezoid rule" begin
        #gaussian PDF should integral to 1.
        pdf(x) = exp(-1/2 * x^2) / sqrt(2π)
        xs = -10:0.1:10
        @test SSSynth.trapezoid_rule(xs, pdf.(xs) * 0.1) - 1.0 < 1e-5
    end
end

@testset "LSF" begin
    wls = 5000:0.35:6000
    R = 1800.0
    flux = zeros(Float64, length(wls))
    flux[500] = 5.0

    convF = SSSynth.constant_R_LSF(flux, wls, R)
    #normalized?
    @test sum(flux) ≈ sum(convF)

    #preserves line center?
    @test argmax(convF) == 500
end
