using Korg, Test, HDF5

@testset "Korg tests" begin

# tools for testing: assert_allclose and assert_allclose_grid
include("utilities.jl") 

# tests for specific parts of the code broken out into their own files. As you add tests, do it 
# this way.
include("cubic_splines.jl")
include("transfer.jl")
include("species.jl")
include("interval.jl")
include("continuum_absorption.jl") # test this after the "Interval" testset
include("partition_funcs.jl")
include("statmech.jl")

@testset "atomic data" begin 
    @test (Korg.MAX_ATOMIC_NUMBER == length(Korg.atomic_masses) == length(Korg.asplund_2009_solar_abundances) 
            == length(Korg.asplund_2020_solar_abundances))
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

@testset "lines" begin
    @testset "linelists" begin 
        @test_throws ArgumentError read_linelist("data/linelists/gfallvac08oct17.stub.dat";
                                                          format="abc")

        @testset "kurucz linelist parsing" begin
            for fname in ["gfallvac08oct17.stub.dat", "gfallvac08oct17-missing-col.stub.dat"]
                kurucz_ll = read_linelist("data/linelists/"*fname, format="kurucz")
                @test issorted(kurucz_ll, by=l->l.wl)
                @test length(kurucz_ll) == 987
                @test kurucz_ll[1].wl ≈ 0.0007234041763337705
                @test kurucz_ll[1].log_gf == -0.826
                @test kurucz_ll[1].species == Korg.species"Be_II"
                @test kurucz_ll[1].E_lower ≈ 17.360339371573698
                @test kurucz_ll[1].gamma_rad ≈ 8.511380382023759e7
                @test kurucz_ll[1].gamma_stark ≈ 0.003890451449942805
                @test kurucz_ll[1].vdW ≈ 1.2302687708123812e-7
            end
        end

        @testset "vald short format, ABO, missing params" begin
            linelist = read_linelist("data/linelists/linelist.vald")
            @test length(linelist) == 6
            @test linelist[1].wl ≈ 3000.0414 * 1e-8
            @test linelist[1].log_gf == -2.957
            @test linelist[1].species == Korg.species"Fe_I"
            @test linelist[1].E_lower ≈ 3.3014
            @test linelist[1].gamma_rad ≈ 1.905460717963248e7
            @test linelist[1].gamma_stark ≈ 0.0001230268770812381
            @test linelist[1].vdW ≈ 4.6773514128719815e-8

            #test imputation of missing broadening parameters
            @test linelist[2].gamma_rad ≈ 818252.5391161365
            @test linelist[2].gamma_stark == linelist[1].gamma_stark
            @test linelist[2].vdW == linelist[1].vdW

            @test linelist[3].gamma_rad == linelist[2].gamma_rad
            @test linelist[3].gamma_stark ≈ 0.00019044182974029873
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
            long_all_cm_air_noquotes = read_linelist("data/linelists/long-extract-all-air-wavenumber-noquotes.vald")
            short_stellar = read_linelist("data/linelists/short-extract-stellar.vald")
            long_stellar = read_linelist("data/linelists/long-extract-stellar.vald")

            @test (length(short_all) == length(short_stellar) == length(long_all_cm_air) == 
                   length(long_all_cm_air_noquotes) == length(long_stellar) == 2)

            #these should be identical since there was no unit conversion
            @test long_all_cm_air[1] == long_all_cm_air_noquotes[1]
            @test short_all[1] == short_stellar[1] == long_stellar[1]

            #when there was unit conversion, they should be approximately equal
            @test short_all[1].wl ≈ long_all_cm_air[1].wl
            @test short_all[1].log_gf == long_all_cm_air[1].log_gf
            @test short_all[1].species == long_all_cm_air[1].species
            @test short_all[1].E_lower ≈ long_all_cm_air[1].E_lower     atol=1e-3
            @test short_all[1].gamma_rad == long_all_cm_air[1].gamma_rad
            @test short_all[1].gamma_stark == long_all_cm_air[1].gamma_stark
            @test short_all[1].vdW == long_all_cm_air[1].vdW
        end

        @testset "vald isotopic scaling" begin
            short_all = read_linelist("data/linelists/isotopic_scaling/short-all-unscaled.vald")
            long_all = read_linelist("data/linelists/isotopic_scaling/long-all-unscaled.vald")
            short_stellar = read_linelist("data/linelists/isotopic_scaling/short-stellar-unscaled.vald")
            long_stellar = read_linelist("data/linelists/isotopic_scaling/long-stellar-unscaled.vald")
            scaled = read_linelist("data/linelists/isotopic_scaling/scaled.vald")
            unscaled_lists = [short_all, long_all, short_stellar, long_stellar]

            for list in unscaled_lists
                for field in [:wl, :species, :E_lower, :gamma_rad, :gamma_stark, :vdW]
                    @test getfield.(scaled, field) == getfield.(list, field)
                end
                for (line1, line2) in zip(scaled, list)
                    @test line1.log_gf .≈ line2.log_gf rtol=1e-2
                end
            end

            #make sure the default isotopic abundances sum to 1.
            for (Z, isotopes) in Korg.isotopic_abundances
                @test sum(values(isotopes)) ≈ 1
            end
        end

        moog_linelist = read_linelist("data/linelists/s5eqw_short.moog"; format="moog")
        @testset "moog linelist parsing" begin
            @test issorted(moog_linelist, by=l->l.wl)
            @test moog_linelist[1].wl ≈ 3729.807 * 1e-8
            @test moog_linelist[1].log_gf ≈ -0.280
            @test moog_linelist[1].species == Korg.species"Ti_I"
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
        wls = (4750 : Δ : 5250) * 1e-8
        Δ *= 1e-8
        amplitude = 7.0
        for σ in [1e-7, 1e-8, 1e-9], γ in [3e-8, 3e-9, 3e-10]
            ϕ = Korg.line_profile.(5e-5, σ, γ, amplitude, wls)
            # the profile isn't perfectly monotonic because the approximation has "seams" at v=5
            # this allows for slight nonmonotonicity
            @test all(diff(ϕ[1:Int(ceil(end/2))]) .> -1e-3*maximum(ϕ))
            @test all(diff(ϕ[Int(ceil(end/2)) : end]) .< 1e-3*maximum(ϕ))
            @test 0.98 < sum(ϕ .* Δ)/amplitude < 1
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

        αs = zeros(length(wls))
        Korg.hydrogen_line_absorption!(αs, wls, 9000.0, 1e11, 1e13, 0.0,
                                       Korg.default_partition_funcs[Korg.species"H_I"](log(9000.0)), 0.0, 15e-7)
        @test assert_allclose_grid(αs_ref, αs, [("λ", wls*1e8, "Å")]; atol=5e-9)

        #make sure that H line absorption doesn't return NaNs on inputs where it used to
        αs = zeros(length(wls))
        wls = 3800 : 0.01 : 4200
        Korg.hydrogen_line_absorption!(αs, wls, 9000.0, 1.1e16, 1, 0.0,
                                       Korg.default_partition_funcs[Korg.species"H_I"](log(9000.0)), 0.0, 15e-7)
        @assert all(.! isnan.(αs))
    end
end

@testset "atmosphere" begin
    @testset "plane-parallel atmosphere" begin
        #the MARCS solar model atmosphere
        atm = Korg.read_model_atmosphere("data/sun.mod")
        @test atm isa Korg.PlanarAtmosphere
        @test length(atm.layers) == 56
        @test issorted([l.temp for l in atm.layers])
        @test atm.layers[1].tau_5000 ≈ 0.00001209483645
        @test atm.layers[1].z == 6.931E+07
        @test atm.layers[1].temp == 4066.8
        @test atm.layers[1].electron_number_density ≈ 3.769664452210607e10
        @test atm.layers[1].number_density ≈ 4.75509171357701e14

        # just make sure these don't error
        Korg.get_tau_5000s(atm)
        Korg.get_zs(atm)
        Korg.get_temps(atm)
        Korg.get_electron_number_densities(atm)
        Korg.get_number_densities(atm)
    end
    @testset "spherical atmosphere" begin
        atm = Korg.read_model_atmosphere(
                "data/s6000_g+1.0_m0.5_t05_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
        @test atm isa Korg.ShellAtmosphere
        @test length(atm.layers) == 56 
        @test issorted([l.temp for l in atm.layers])
        @test atm.R == 2.5827E+12
        @test atm.layers[1].tau_5000 ≈ 4.584584692493259e-5
        @test atm.layers[1].z == 2.222e11
        @test atm.layers[1].temp == 3935.2
        @test atm.layers[1].electron_number_density ≈ 1.7336231777439526e8
        @test atm.layers[1].number_density ≈ 1.5411190391302566e12
    end

    @testset "atmosphere type conversion" begin
        atm = Korg.read_model_atmosphere("data/sun.mod")
        atm2 = Korg.PlanarAtmosphere(Korg.ShellAtmosphere(atm, 7e10)) #arbitrary radius
        @test atm.layers == atm2.layers

        atm = Korg.read_model_atmosphere(
                "data/s6000_g+1.0_m0.5_t05_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
        atm2 = Korg.ShellAtmosphere(Korg.PlanarAtmosphere(atm), 1.0)
        @test [l.tau_5000 for l in atm.layers]                == [l.tau_5000 for l in atm2.layers]
        @test [l.z for l in atm.layers]                       == [l.z for l in atm2.layers]
        @test [l.temp for l in atm.layers]                    == [l.temp for l in atm2.layers]
        @test [l.number_density for l in atm.layers]          == [l.number_density for l in atm2.layers]
        @test [l.electron_number_density for l in atm.layers] == [l.electron_number_density for l in atm2.layers]
    end
end

@testset "synthesis" begin

    @testset "abundances" begin
        @test_throws ArgumentError format_A_X(0.0, Dict("H"=>13))
        @test (format_A_X() 
                == format_A_X(0)
                == format_A_X(Dict{String, Float64}())
                == format_A_X(Dict{Int, Float64}()))

        @testset for metallicity in [0.0, 0.5], abundances in [Dict(), Dict("C"=>1.1)], solar_relative in [true, false]
            A_X = format_A_X(metallicity, abundances; 
                                       solar_abundances=Korg.asplund_2020_solar_abundances,
                                       solar_relative=solar_relative)

            #correct absolute abundances?
            if "C" in keys(abundances)
                if solar_relative
                    @test A_X[6] ≈ Korg.asplund_2020_solar_abundances[6] + 1.1
                else
                    @test A_X[6] ≈ 1.1
                end
            end
            @test A_X[7:end] ≈ Korg.asplund_2020_solar_abundances[7:end] .+ metallicity
            @test A_X[1:2] == Korg.asplund_2020_solar_abundances[1:2]
        end
    end
end

@testset "LSF" begin
    wls = 5000:0.35:6000
    R = 1800.0
    flux = zeros(Float64, length(wls))
    flux[500] = 5.0

    convF = Korg.constant_R_LSF(flux, wls, R)
    convF_4sigma = Korg.constant_R_LSF(flux, wls, R; window_size=4)

    #normalized?
    @test sum(flux) ≈ sum(convF)
    @test sum(flux) ≈ sum(convF_4sigma)

    #preserves line center?
    @test argmax(convF) == 500
    @test argmax(convF_4sigma) == 500

    #make sure the window_size argument is doing something
    @test !(convF ≈ convF_4sigma)
end

@testset "air <--> vacuum" begin
    wls = collect(2000.0:π:10000.0)
    @test Korg.vacuum_to_air.(Korg.air_to_vacuum.(wls)) ≈ wls rtol=1e-3
    @test Korg.air_to_vacuum.(Korg.vacuum_to_air.(wls)) ≈ wls rtol=1e-3

    #try it in cgs
    wls .*= 1e8
    @test Korg.vacuum_to_air.(Korg.air_to_vacuum.(wls)) ≈ wls rtol=1e-3
    @test Korg.air_to_vacuum.(Korg.vacuum_to_air.(wls)) ≈ wls rtol=1e-3

    #units should be automatically chosen
    @test Korg.vacuum_to_air.(Korg.air_to_vacuum.(wls*1e-8)*1e8) ≈ wls rtol=1e-3
    @test Korg.air_to_vacuum.(Korg.vacuum_to_air.(wls*1e-8)*1e8) ≈ wls rtol=1e-3
    @test Korg.vacuum_to_air.(Korg.air_to_vacuum.(wls)*1e8)*1e-8 ≈ wls rtol=1e-3
    @test Korg.air_to_vacuum.(Korg.vacuum_to_air.(wls)*1e8)*1e-8 ≈ wls rtol=1e-3
end

@testset "synthesize wavelength handling" begin
    atm = read_model_atmosphere("data/sun.mod")
    wls = 15000:0.01:15500
    A_X = format_A_X()
    @test synthesize(atm, [], A_X, 15000, 15500).wavelengths ≈ wls
    @test synthesize(atm, [], A_X, 15000, 15500; air_wavelengths=true).wavelengths ≈ Korg.air_to_vacuum.(wls)
    @test_throws ArgumentError synthesize(atm, [], A_X, 15000, 15500; air_wavelengths=true, 
                                          wavelength_conversion_warn_threshold=1e-20)
    @test_throws ArgumentError synthesize(atm, [], A_X, 2000, 8000, air_wavelengths=true)
end

@testset "line buffer" begin
    #strong line at 4999 Å
    line1 = Korg.Line(4999e-8, 1.0, Korg.species"Na I", 0.0)
    #strong line at 4997 Å
    line2 = Korg.Line(4997e-8, 1.0, Korg.species"Na I", 0.0)
    atm = read_model_atmosphere("data/sun.mod")

    #use a 2 Å line buffer so only line1 in included
    sol_no_lines = synthesize(atm, [], format_A_X(), 5000, 5000; line_buffer=2.0) #synthesize at 5000 Å only
    sol_one_lines = synthesize(atm, [line1], format_A_X(), 5000, 5000; line_buffer=2.0) 
    sol_two_lines = synthesize(atm, [line1, line2], format_A_X(), 5000, 5000; line_buffer=2.0) 

    @test sol_no_lines.flux != sol_one_lines.flux
    @test sol_two_lines.flux == sol_one_lines.flux
end

@testset "autodiff" begin
    using ForwardDiff

    linelist = read_linelist("data/linelists/5000-5005.vald")
    wls = 6564:0.01:6565
    for atm_file in ["data/sun.mod",
             "data/s6000_g+1.0_m0.5_t05_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod"]
        atm = read_model_atmosphere(atm_file)
        flux(p) = synthesize(atm, linelist, format_A_X(p[1], Dict("Ni"=>p[2])), 
                             wls; vmic=p[3]).flux
        #make sure this works.
        J = ForwardDiff.jacobian(flux, [0.0, 0.0, 1.5])
        @test .! any(isnan.(J))
    end
end

end #top-level testset