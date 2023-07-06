using Korg, Test, HDF5, ForwardDiff, FiniteDiff

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
include("linelist.jl")
include("fit.jl")

@testset "atomic data" begin 
    @test (Korg.MAX_ATOMIC_NUMBER 
            == length(Korg.atomic_masses) 
            == length(Korg.asplund_2009_solar_abundances) 
            == length(Korg.asplund_2020_solar_abundances) 
            == length(Korg.grevesse_2007_solar_abundances) 
            == length(Korg.magg_2022_solar_abundances))

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

    a = 1:10
    @test Korg.move_bounds(a, 0, 0, 5.5, 0.1) == (6, 5)

    a = [3:5, 11:0.5:12.5, 16:20]
    @test Korg.move_bounds(a, 0, 0, -1, 1) == (1, 0)
    @test Korg.move_bounds(a, 0, 0, 3, 1) == (1, 2)
    @test Korg.move_bounds(a, 0, 0, 5, 6) == (1, 4)
    @test Korg.move_bounds(a, 0, 0, 12.5, 0.6) == (6, 7)
    @test Korg.move_bounds(a, 0, 0, 50, 5) == (1, 0)
end

@testset "line profiles" begin
    @testset "generic line profile" begin
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

        fid = h5open(fname) 
        T = HDF5.read_attribute(fid["profile"], "T")
        ne = HDF5.read_attribute(fid["profile"], "ne")
       nH_I = HDF5.read_attribute(fid["profile"], "nH_I")
        wls = (HDF5.read_attribute(fid["profile"], "start_wl") :
               HDF5.read_attribute(fid["profile"], "wl_step") : 
               HDF5.read_attribute(fid["profile"], "stop_wl") )
        close(fid)

        αs = zeros(length(wls))
        Korg.hydrogen_line_absorption!(αs, [wls], 9000.0, ne, nH_I, 0.0,
                                       Korg.default_partition_funcs[Korg.species"H_I"](log(9000.0)), 
                                       0.0, 15e-7, use_MHD=false) 
        @test assert_allclose_grid(αs_ref, αs, [("λ", wls*1e8, "Å")]; atol=5e-9)

        #make sure that H line absorption doesn't return NaNs on inputs where it used to
        αs = zeros(length(wls))
        wls = 3800 : 0.01 : 4200
        Korg.hydrogen_line_absorption!(αs, [wls], 9000.0, 1.1e16, 1, 0.0,
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
        @test (format_A_X() 
                == format_A_X(0)
                == format_A_X(0, 0)
                == format_A_X(Dict{String, Float64}())
                == format_A_X(Dict{Int, Float64}())
                == format_A_X(0, Dict(1=>0.0); solar_relative=true)
                == format_A_X(0, 0, Dict(1=>0.0); solar_relative=true)
                == format_A_X(0, Dict("H"=>0.0); solar_relative=true)
                == format_A_X(0, Dict(1=>12.0); solar_relative=false)
                == format_A_X(0, Dict("H"=>12.0); solar_relative=false))
        
        # make sure silly H abundances are caught
        @test_throws ArgumentError format_A_X(0.0, Dict("H"=>0); solar_relative=false)
        @test_throws ArgumentError format_A_X(0.0, Dict(1=>0); solar_relative=false)
        @test_throws ArgumentError format_A_X(0.0, Dict("H"=>12); solar_relative=true)
        @test_throws ArgumentError format_A_X(0.0, Dict(1=>12); solar_relative=true)

        @test Korg.get_alpha_H(format_A_X(0.1)) ≈ 0.1 atol=1e-6
        @test Korg.get_alpha_H(format_A_X(0.0, 0.1)) ≈ 0.1 atol=1e-6
        @test Korg.get_alpha_H(format_A_X(-0.2)) ≈ -0.2 atol=1e-6
        @test Korg.get_alpha_H(format_A_X(-2, -0.2)) ≈ -0.2 atol=1e-6
        @test Korg.get_metals_H(format_A_X(0.1)) ≈ 0.1 atol=1e-6
        @test Korg.get_metals_H(format_A_X(-0.2)) ≈ -0.2 atol=1e-6
        @test Korg.get_metals_H(Korg.grevesse_2007_solar_abundances; 
                                solar_abundances=Korg.grevesse_2007_solar_abundances) ≈ 0 atol=1e-6
        @test Korg.get_alpha_H(Korg.grevesse_2007_solar_abundances;
                               solar_abundances=Korg.grevesse_2007_solar_abundances) ≈ 0 atol=1e-6

        @test format_A_X(1.1) != format_A_X(1.1, 0)
        @test format_A_X(1.1)[50] == format_A_X(1.1, 0)[50] == format_A_X(-1, -2, Dict(50=>1.1))[50]

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
    wls = 5900:0.35:6100
    R = 1800.0
    spike_ind = 286 #in the middle
    flux = zeros(Float64, length(wls))
    flux[spike_ind] = 5.0

    convF = Korg.constant_R_LSF(flux, wls, R)
    convF_5sigma = Korg.constant_R_LSF(flux, wls, R; window_size=5)
    convF_mat = Korg.compute_LSF_matrix(wls, wls, R) * flux
    convF_mat5 = Korg.compute_LSF_matrix(wls, wls, R; window_size=5) * flux
    downsampled_wls = 5950:0.4:6050
    convF_mat_downsample = Korg.compute_LSF_matrix(wls, downsampled_wls, R; window_size=5) * flux

    # normalized?
    @test sum(flux) ≈ sum(convF)  rtol=1e-3
    @test sum(flux) ≈ sum(convF_5sigma) rtol=1e-3
    @test sum(flux) ≈ sum(convF_mat) rtol=1e-3
    @test sum(flux) ≈ sum(convF_mat5) rtol=1e-3
    @test sum(flux) ≈ sum(convF_mat5) rtol=1e-3
    @test sum(flux) * step(wls) ≈ sum(convF_mat_downsample) * step(downsampled_wls) rtol=1e-3
    

    # preserves line center?
    @test argmax(convF) == spike_ind
    @test argmax(convF_5sigma) == spike_ind
    @test argmax(convF_mat) == spike_ind
    @test argmax(convF_mat5) == spike_ind

    # make sure the default window_size values are OK
    @test assert_allclose(convF, convF_mat5; atol=1e-4)
    @test assert_allclose(convF_5sigma, convF_mat5; atol=1e-4)
    @test assert_allclose(convF_mat, convF_mat5; atol=1e-4)

    # but also check that they are definitely doing something
    @test !(convF ≈ convF_5sigma)
    @test !(convF_mat ≈ convF_mat5)
end

@testset "rotation" begin
    wls = 4090:0.01:5010
    flux = zeros(length(wls))
    flux[990:1010] .= 1

    for vsini in [0.0, 1.0, 5.0, 10.0, 20.0], ε in [0.1, 0.6, 0.9]
        # also test handling of multiple wl ranges
        for wls in [wls, [4090:0.01:5007, 5007.01:0.01:5010]]
            rflux = Korg.apply_rotation(flux, wls, vsini, ε)
            rflux2 = Korg.apply_rotation(flux, wls * 1e-8, vsini, ε)

            # rotational kernel is normalized
            @test sum(flux) ≈ sum(rflux) rtol=1e-2
            @test sum(flux) ≈ sum(rflux2) rtol=1e-2

            @test rflux == rflux2 # wl units shouldn't matter
        end
    end
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

    # test multiple line windows
    r1 = 5000:0.01:5001
    r2 = 6000:0.01:6001
    sol1 = synthesize(atm, [], A_X, [r1]; hydrogen_lines=true)
    sol2 = synthesize(atm, [], A_X, [r2]; hydrogen_lines=true)
    sol3 = synthesize(atm, [], A_X, [r1, r2]; hydrogen_lines=true)

    @test sol1.wavelengths == sol3.wavelengths[sol3.subspectra[1]]
    @test sol2.wavelengths == sol3.wavelengths[sol3.subspectra[2]]
    @test sol1.flux == sol3.flux[sol3.subspectra[1]]
    @test sol2.flux == sol3.flux[sol3.subspectra[2]]
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

@testset "continuum" begin
    atm = read_model_atmosphere("data/sun.mod") 
    linelist = [Korg.Line(6560, 1.0, Korg.species"Na I", 0.0)]
    A_X = format_A_X()
    sol = synthesize(atm, linelist, A_X, 6555, 6560)
    cntm_sol = synthesize(atm, [], A_X, 6555, 6560; hydrogen_lines=false, return_cntm=false)
    @test all(sol.cntm .≈ cntm_sol.flux)
    @test isnothing(sol.cntm)
end

@testset "autodiff" begin
    using ForwardDiff, FiniteDiff

    linelist = read_linelist("data/linelists/5000-5005.vald")
    wls = [6564:0.01:6565]
    for (atm_file, threshold) in [("data/sun.mod", 0.1),
             ("data/s6000_g+1.0_m0.5_t05_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod", 100)]
        atm = read_model_atmosphere(atm_file)
        # the second model atmosphere happens to be in a weird (probably unphysical) part of 
        # parameter space where the electron number densities calculated doesn't match the marcs
        # numbers.  We suppress the warnings in this case by setting the threshold to 100.
        flux(p) = synthesize(atm, linelist, format_A_X(p[1], Dict("Ni"=>p[2])), wls; 
                             vmic=p[3], electron_number_density_warn_threshold=threshold).flux
        #make sure this works.
        J = ForwardDiff.jacobian(flux, [0.0, 0.0, 1.5])
        @test .! any(isnan.(J))
    end

    @testset "line params" begin
        atm = read_model_atmosphere("data/sun.mod")
        function f(loggf)
            linelist = [Korg.Line(5000e-8, loggf, Korg.species"Na I", 0.0)]
            synthesize(atm, linelist, format_A_X(), 5000, 5000).flux[1]
        end
        @test FiniteDiff.finite_difference_derivative(f, 0.0) ≈ ForwardDiff.derivative(f, 0.0) rtol=1e-4
    end

end

end #top-level testset
