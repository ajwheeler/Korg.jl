using Korg, Test, HDF5

include("continuum_opacity.jl")

@testset "atomic data" begin 
    @test Korg.Natoms == length(Korg.atomic_masses) == length(Korg.solar_abundances)
    @test (Korg.get_mass(Korg.Formula("CO")) ≈ 
           Korg.get_mass(Korg.Formula("C")) + Korg.get_mass(Korg.Formula("O")))
    @test Korg.get_mass(Korg.Formula("C2")) ≈ 2Korg.get_mass(Korg.Formula("C"))
end

@testset "ionization energies" begin
    @test length(Korg.ionization_energies) == 92
    @test Korg.ionization_energies[Korg.atomic_numbers["H"]] == [13.5984, -1.000, -1.000]
    @test Korg.ionization_energies[Korg.atomic_numbers["Ru"]] == [7.3605, 16.760, 28.470]
    @test Korg.ionization_energies[Korg.atomic_numbers["U"]] == [6.1940, 11.590, 19.800]
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
    coef = 2.0 * (2.0*π*Korg.electron_mass_cgs*Korg.kboltz_cgs / Korg.hplanck_cgs^2)^1.5
    RHS = coef * T^1.5 * exp(-Korg.RydbergH_eV/(Korg.kboltz_eV*T))/H_I_partition_val
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
    @test issorted(Korg.partition_funcs[Korg.Species("O_I")].(Ts))
    @test issorted(Korg.partition_funcs[Korg.Species("O_II")].(Ts))
    @test issorted(Korg.partition_funcs[Korg.Species("O_III")].(Ts))
    @test issorted(Korg.partition_funcs[Korg.Species("CN_I")].(Ts))
end

@testset "stat mech" begin
    @testset "pure Hydrogen atmosphere" begin
        nH_tot = 1e15
        # specify χs and Us to decouple this testset from other parts of the code
        χs = Dict(1=>[Korg.RydbergH_eV, -1.0, -1.0])
        Us = Dict([Korg.Species("H_I")=>(T -> 2.0), Korg.Species("H_II")=>(T -> 1.0)])
        # iterate from less than 1% ionized to more than 99% ionized
        for T in [3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 1e4, 1.1e4, 1.2e4, 1.3e4, 1.4e4, 1.5e5]
            nₑ = electron_ndens_Hplasma(nH_tot, T, 2.0)
            wII, wIII = Korg.saha_ion_weights(T, nₑ, 1, χs, Us)
            @test wIII ≈ 0.0 rtol = 1e-15
            rtol = (T == 1.5e5) ? 1e-9 : 1e-14
            @test wII/(1 + wII + wIII) ≈ (nₑ/nH_tot) rtol= rtol
        end
    end

    @testset "monotonic N ions Temperature dependence" begin
        weights = [Korg.saha_ion_weights(T, 1.0, 7, Korg.ionization_energies, 
                                            Korg.partition_funcs) for T in 1:100:10000]
        #N II + NIII grows with T === N I shrinks with T
        @test issorted(first.(weights) + last.(weights))
        
        # NIII grows with T
        @test issorted(last.(weights))
    end

    @testset "molecular equilibrium" begin
        #solar abundances
        abundances = Korg.get_absolute_abundances(0.0, Dict())
        nₜ = 1e15 
        nₑ = 1e-3 * nₜ #arbitrary

        MEQs = Korg.molecular_equilibrium_equations(abundances, Korg.ionization_energies, 
                                                       Korg.partition_funcs, 
                                                       Korg.equilibrium_constants)
        @test MEQs.atoms == 0x01:Korg.Natoms

        n = Korg.molecular_equilibrium(MEQs, 5700.0, nₜ, nₑ)
        #make sure number densities are sensible
        @test (n[Korg.Species("C_III")] < n[Korg.Species("C_II")] < n[Korg.Species("C_I")] < 
               n[Korg.Species("H_II")] < n[Korg.Species("H_I")])

        #total number of carbons is correct
        total_C = map(collect(keys(n))) do species
            if species == Korg.Species("C2")
                n[species] * 2
            elseif 0x06 in Korg.get_atoms(species.formula)
                n[species]
            else
                0.0
            end
        end |> sum
        @test total_C ≈ abundances[Korg.atomic_numbers["C"]] * nₜ
    end
end

@testset "lines" begin
    @testset "linelists" begin 
        @testset "species codes" begin
            @test Korg.Species("01.00")   == Korg.Species("H I")
            @test Korg.Species("101.0")   == Korg.Species("H2 I")
            @test Korg.Species("01.0000") == Korg.Species("H I")
            @test Korg.Species("02.01")   == Korg.Species("He II")
            @test Korg.Species("02.1000") == Korg.Species("He II")
            @test Korg.Species("0608")    == Korg.Species("CO I")
            @test Korg.Species("0606")    == Korg.Species("C2 I")
            @test Korg.Species("606")     == Korg.Species("C2 I")
            @test Korg.Species("0608.00") == Korg.Species("CO I")
            @test Korg.Species("OOO")     == Korg.Species("O3")

            @test_throws ArgumentError Korg.Species("06.05.04")
            @test_throws Exception Korg.Species("99.01")
        end

        @testset "distinguish atoms from molecules" begin
            @test Korg.ismolecule(Korg.Formula("H2"))
            @test Korg.ismolecule(Korg.Formula("CO"))
            @test !Korg.ismolecule(Korg.Formula("H"))
            @test !Korg.ismolecule(Korg.Formula("Li"))
        end

        @testset "break molecules into atoms" begin
            @test Korg.get_atoms(Korg.Formula("CO")) == [0x06, 0x08]
            @test Korg.get_atoms(Korg.Formula("C2")) == [0x06, 0x06]
            @test Korg.get_atoms(Korg.Formula("MgO")) == [0x08, 0x0c]
        end

        @test_throws ArgumentError read_linelist("data/linelists/gfallvac08oct17.stub.dat";
                                                          format="abc")
        @test_throws ArgumentError read_linelist("data/linelists/no-isotopic-scaling.vald")

        kurucz_ll = read_linelist("data/linelists/gfallvac08oct17.stub.dat", format="kurucz")
        @testset "kurucz linelist parsing" begin
            @test issorted(kurucz_ll, by=l->l.wl)
            @test length(kurucz_ll) == 987
            @test kurucz_ll[1].wl ≈ 72320.699 * 1e-8
            @test kurucz_ll[1].log_gf == -0.826
            @test kurucz_ll[1].species == Korg.Species("Be_II")
            @test kurucz_ll[1].E_lower ≈ 17.360339371573698
            @test kurucz_ll[1].gamma_rad ≈ 8.511380382023759e7
            @test kurucz_ll[1].gamma_stark ≈ 0.003890451449942805
            @test kurucz_ll[1].vdW ≈ 1.2302687708123812e-7
        end

        @testset "vald short format, ABO, missing params" begin
            linelist = read_linelist("data/linelists/linelist.vald")
            @test length(linelist) == 6
            @test linelist[1].wl ≈ 3000.0414 * 1e-8
            @test linelist[1].log_gf == -2.957
            @test linelist[1].species == Korg.Species("Fe_I")
            @test linelist[1].E_lower ≈ 3.3014
            @test linelist[1].gamma_rad ≈ 1.905460717963248e7
            @test linelist[1].gamma_stark ≈ 0.0001230268770812381
            @test linelist[1].vdW ≈ 4.6773514128719815e-8

            #test imputation of missing broadening parameters
            @test linelist[2].gamma_rad ≈ 818252.5391161365
            @test linelist[2].gamma_stark == linelist[1].gamma_stark
            @test linelist[2].vdW == linelist[1].vdW

            @test linelist[3].gamma_rad == linelist[2].gamma_rad
            @test linelist[3].gamma_stark == 5.848503287015111e-25
            @test linelist[3].vdW == linelist[1].vdW

            @test linelist[4].gamma_rad == linelist[1].gamma_rad
            @test linelist[4].gamma_stark == linelist[3].gamma_stark
            @test linelist[4].vdW == 9.953360714197118e-8

            @test linelist[5].gamma_rad == linelist[1].gamma_rad
            @test linelist[5].gamma_stark == linelist[1].gamma_stark
            @test linelist[5].vdW == linelist[4].vdW

            #ABO params
            @test linelist[6].vdW[1] ≈ 1.3917417470792187e-14
            @test linelist[6].vdW[2] ≈ 0.227
        end

        @testset "vald various formats" begin
            short_all = read_linelist("data/linelists/short-extract-all.vald")
            long_all_cm_air = read_linelist("data/linelists/long-extract-all-air-wavenumber.vald")
            short_stellar = read_linelist("data/linelists/short-extract-stellar.vald")
            long_stellar = read_linelist("data/linelists/long-extract-stellar.vald")

            @test (length(short_all) == length(short_stellar) == length(long_all_cm_air) == 
                   length(long_stellar) == 1)

            @test short_all[1] == short_stellar[1] == long_stellar[1]

            @test short_all[1].wl ≈ long_all_cm_air[1].wl
            @test short_all[1].log_gf == long_all_cm_air[1].log_gf
            @test short_all[1].species == long_all_cm_air[1].species
            @test short_all[1].E_lower ≈ long_all_cm_air[1].E_lower     atol=1e-3
            @test short_all[1].gamma_rad == long_all_cm_air[1].gamma_rad
            @test short_all[1].gamma_stark == long_all_cm_air[1].gamma_stark
            @test short_all[1].vdW == long_all_cm_air[1].vdW
        end

        moog_linelist = read_linelist("data/linelists/s5eqw_short.moog"; format="moog")
        @testset "moog linelist parsing" begin
            @test issorted(moog_linelist, by=l->l.wl)
            @test moog_linelist[1].wl ≈ 3729.807 * 1e-8
            @test moog_linelist[1].log_gf ≈ -0.280
            @test moog_linelist[1].species == Korg.Species("Ti_I")
            @test moog_linelist[2].E_lower ≈ 3.265
        end
    end

    @testset "move_bounds" begin
        a = 0.5 .+ (1:9)
        for lb in [1, 3, 9], ub in [1, 5, 9]
            @test Korg.move_bounds(a, lb, ub, 5., 2.) == (3, 6)
            @test Korg.move_bounds(a, lb, ub, 0., 3.) == (1, 2)
            @test Korg.move_bounds(a, lb, ub, 6., 4.) == (2, 9)
            @test Korg.move_bounds(collect(a), lb, ub, 5., 2.) == (3, 6)
            @test Korg.move_bounds(collect(a), lb, ub, 0., 3.) == (1, 2)
            @test Korg.move_bounds(collect(a), lb, ub, 6., 4.) == (2, 9)
        end
    end

    @testset "line profile" begin
        Δ = 0.01
        wls = (4955 : Δ : 5045) * 1e-8
        Δ *= 1e-8
        amplitude = 7.0
        for Δλ_D in [1e-7, 1e-8, 1e-9], Δλ_L in [1e-8, 1e-9]
            ϕ = Korg.line_profile.(5e-5, 1/Δλ_D, Δλ_L, amplitude, wls)
            @test issorted(ϕ[1 : Int(ceil(end/2))])
            @test issorted(ϕ[Int(ceil(end/2)) : end], rev=true)
            @test 0.99 < sum(ϕ .* Δ)/amplitude < 1
        end
    end

    @testset "hydrogen stark profiles" begin
        # This test data was generated with Korg.hydrogen_line_absorption shortly
        # after writing the function. This data is consistent with the results
        # produced by the Fortran code distributed with Stehle & Hutcheon 1999
        fname = "data/lyman_absorption.h5"
        αs_ref = h5read(fname,  "profile")

        fid = h5open("data/lyman_absorption.h5") 
        T = HDF5.read_attribute(fid["profile"], "T")
        ne = HDF5.read_attribute(fid["profile"], "ne")
        nH_I = HDF5.read_attribute(fid["profile"], "nH_I")
        wls = (HDF5.read_attribute(fid["profile"], "start_wl") :
               HDF5.read_attribute(fid["profile"], "wl_step") : 
               HDF5.read_attribute(fid["profile"], "stop_wl") )
        close(fid)

        αs = Korg.hydrogen_line_absorption(wls, 9000.0, 1e11, 1e13, 
                                           Korg.partition_funcs[Korg.Species("H_I")], 
                                           Korg.hline_stark_profiles, 0.0)
        @test αs_ref ≈ αs rtol=1e-5
    end
end

@testset "atmosphere" begin
    #the MARCS solar model atmosphere
    atmosphere = Korg.read_model_atmosphere("data/sun.krz")
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
        @test_throws ArgumentError Korg.get_absolute_abundances(0.0, Dict("H"=>13))

        @testset for metallicity in [0.0, 1.0], A_X in [Dict(), Dict("C"=>9)]
            nxnt = Korg.get_absolute_abundances(metallicity, A_X)

            #correct absolute abundances?
            if "C" in keys(A_X)
                @test log10(nxnt[6]/nxnt[1]) + 12 ≈ 9
            end
            @test log10(nxnt[2]/nxnt[1]) + 12 ≈ Korg.solar_abundances[2]
            @test log10(nxnt[Korg.atomic_numbers["Ba"]]/nxnt[1]) + 12 ≈ 
                Korg.solar_abundances[Korg.atomic_numbers["Ba"]] + metallicity

            #normalized?
            @test sum(nxnt) ≈ 1
        end
    end

    @testset "trapezoid rule" begin
        #gaussian PDF should integral to 1.
        pdf(x) = exp(-1/2 * x^2) / sqrt(2π)
        xs = -10:0.1:10
        @test Korg.trapezoid_rule(xs, pdf.(xs) * 0.1) - 1.0 < 1e-5
    end
end

@testset "LSF" begin
    wls = 5000:0.35:6000
    R = 1800.0
    flux = zeros(Float64, length(wls))
    flux[500] = 5.0

    convF = Korg.constant_R_LSF(flux, wls, R)
    #normalized?
    @test sum(flux) ≈ sum(convF)

    #preserves line center?
    @test argmax(convF) == 500
end

@testset "air <--> vacuum" begin
    wls = collect(2000.0:π:10000.0)
    @test vacuum_to_air.(air_to_vacuum.(wls)) ≈ wls rtol=1e-3
    @test air_to_vacuum.(vacuum_to_air.(wls)) ≈ wls rtol=1e-3

    #try it in cgs
    wls .*= 1e8
    @test vacuum_to_air.(air_to_vacuum.(wls)) ≈ wls rtol=1e-3
    @test air_to_vacuum.(vacuum_to_air.(wls)) ≈ wls rtol=1e-3

    #units should be automatically chosen
    @test vacuum_to_air.(air_to_vacuum.(wls*1e-8)*1e8) ≈ wls rtol=1e-3
    @test air_to_vacuum.(vacuum_to_air.(wls*1e-8)*1e8) ≈ wls rtol=1e-3
    @test vacuum_to_air.(air_to_vacuum.(wls)*1e8)*1e-8 ≈ wls rtol=1e-3
    @test air_to_vacuum.(vacuum_to_air.(wls)*1e8)*1e-8 ≈ wls rtol=1e-3
end

@testset "autodiff" begin
    using ForwardDiff

    atm = read_model_atmosphere("data/sun.krz")
    linelist = read_linelist("data/linelists/5000-5005.vald")
    wls = 5000:0.01:5005
    flux(p) = synthesize(atm, linelist, wls; metallicity=p[1], abundances=Dict(["Ni"=>p[2]]), 
                         vmic=p[3]).flux

    #make sure this works.
    ∇f = ForwardDiff.jacobian(flux, [0.0, 0.0, 1.5])
end
