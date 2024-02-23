using Korg, Test, HDF5, ForwardDiff, FiniteDiff

@testset "Korg tests" begin

# tools for testing: assert_allclose and assert_allclose_grid
include("utilities.jl") 

# tests for specific parts of the code broken out into their own files. As you add tests, do it 
# this way.
include("molecular_cross_sections.jl")
include("cubic_splines.jl")
include("transfer.jl")
include("species.jl")
include("interval.jl")
include("continuum_absorption.jl") # test this after the "Interval" testset
include("partition_funcs.jl")
include("statmech.jl")
include("linelist.jl")
include("fit.jl")
include("autodiff.jl")
include("autodiffable_conv.jl")
include("atmosphere.jl")

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
    @test Korg.move_bounds(a, 0, 0, 50, 5) == (13, 12)

    # check that indices are appropriately out of order when the range is outside a
    for a in [1:10, collect(1:10), [1:5, 6:10]]
        @test Korg.move_bounds(a, 1, 1, 13.5, 0.2) == (11, 10)
        @test Korg.move_bounds(a, 1, 1, -5, 0.2) == (1, 0)
    end
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
        Korg.hydrogen_line_absorption!(αs, wls, [wls], 9000.0, ne, nH_I, 0.0,
                                       Korg.default_partition_funcs[Korg.species"H_I"](log(9000.0)), 
                                       0.0, 15e-7, use_MHD=false) 
        @test assert_allclose_grid(αs_ref, αs, [("λ", wls*1e8, "Å")]; atol=5e-9)

        #make sure that H line absorption doesn't return NaNs on inputs where it used to
        αs = zeros(length(wls))
        wls = 3800 : 0.01 : 4200
        Korg.hydrogen_line_absorption!(αs, wls, [wls], 9000.0, 1.1e16, 1, 0.0,
                                       Korg.default_partition_funcs[Korg.species"H_I"](log(9000.0)), 0.0, 15e-7)
        @assert all(.! isnan.(αs))
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

        atol = 1e-5
        @test Korg.get_alpha_H(format_A_X(0.1)) ≈ 0.1 atol=atol
        @test Korg.get_alpha_H(format_A_X(0.0, 0.1)) ≈ 0.1 atol=atol
        @test Korg.get_alpha_H(format_A_X(-0.2)) ≈ -0.2 atol=atol
        @test Korg.get_alpha_H(format_A_X(-2, -0.2)) ≈ -0.2 atol=atol
        @test Korg.get_metals_H(format_A_X(0.1)) ≈ 0.1 atol=atol
        @test Korg.get_metals_H(format_A_X(-0.2)) ≈ -0.2 atol=atol
        @test Korg.get_metals_H(format_A_X(0.1, 0.5)) ≈ 0.1 atol=atol
        @test Korg.get_metals_H(format_A_X(-0.2, 0.5)) ≈ -0.2 atol=atol
        @test Korg.get_metals_H(Korg.grevesse_2007_solar_abundances; 
                                solar_abundances=Korg.grevesse_2007_solar_abundances) ≈ 0 atol=atol
        @test Korg.get_alpha_H(Korg.grevesse_2007_solar_abundances;
                               solar_abundances=Korg.grevesse_2007_solar_abundances) ≈ 0 atol=atol

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
    flux = zeros(length(wls))
    flux[spike_ind] = 5.0

    # should't matter if you split wls into 2 ranges
    @test Korg.compute_LSF_matrix([5900:0.35:5935, 5935.35:0.35:6100], wls, R) ≈ Korg.compute_LSF_matrix(wls, wls, R) 

    # should't matter if R is a function
    @test Korg.compute_LSF_matrix(wls, wls, wl->R) ≈ Korg.compute_LSF_matrix(wls, wls, R) 

    # synth_wls can be a slightly fuzzy grid
    obs_wls = 5000 : 1.0 : 5010
    range_wls = 5000 : 0.01 : 5010
    array_wls = collect(range_wls) # range -> array
    array_wls[2:end-1] .+= 1e-8 * ones(length(range_wls)-2) # slighly perturb (except the endpoints)
    @test Korg.compute_LSF_matrix(array_wls, obs_wls, R) ≈ Korg.compute_LSF_matrix(range_wls, obs_wls, R) 
    @test_throws ArgumentError Korg.compute_LSF_matrix(array_wls, obs_wls, R; step_tolerance=1e-9)

    convF = Korg.apply_LSF(flux, wls, R)
    convF_5sigma = Korg.apply_LSF(flux, wls, R; window_size=5)
    convF_mat = Korg.compute_LSF_matrix(wls, wls, R) * flux
    convF_mat5 = Korg.compute_LSF_matrix(wls, wls, R; window_size=5) * flux
    convF_mat_vec = Korg.compute_LSF_matrix([5900:0.35:5935, 5935.35:0.35:6100], wls, R) * flux
    convF_changing_R = Korg.apply_LSF(flux, wls, wl->wl/6000*R)
    convF_mat_changing_R = Korg.compute_LSF_matrix(wls, wls, wl->wl/6000*R) * flux

    @test convF_mat_changing_R ≈ convF_changing_R rtol=1e-10

    downsampled_wls = 5950:0.4:6050
    convF_mat_downsample = Korg.compute_LSF_matrix(wls, downsampled_wls, R; window_size=5) * flux

    # normalized?
    @test sum(flux) ≈ sum(convF)  rtol=1e-3
    @test sum(flux) ≈ sum(convF_5sigma) rtol=1e-3
    @test sum(flux) ≈ sum(convF_mat) rtol=1e-3
    @test sum(flux) ≈ sum(convF_mat5) rtol=1e-3
    @test sum(flux) ≈ sum(convF_mat5) rtol=1e-3
    @test sum(flux) ≈ sum(convF_mat_vec) rtol=1e-3
    @test sum(flux) ≈ sum(convF_changing_R) rtol=1e-3
    @test sum(flux) ≈ sum(convF_mat_changing_R) rtol=1e-3
    @test sum(flux) * step(wls) ≈ sum(convF_mat_downsample) * step(downsampled_wls) rtol=1e-3

    # preserves line center?
    @test argmax(convF) == spike_ind
    @test argmax(convF_5sigma) == spike_ind
    @test argmax(convF_mat) == spike_ind
    @test argmax(convF_mat5) == spike_ind
    @test argmax(convF_mat_vec) == spike_ind
    @test argmax(convF_changing_R) == spike_ind
    @test argmax(convF_mat_changing_R) == spike_ind

    # make sure the default window_size values are OK
    @test assert_allclose(convF, convF_mat5; atol=1e-4)
    @test assert_allclose(convF_5sigma, convF_mat5; atol=1e-4)
    @test assert_allclose(convF_mat, convF_mat5; atol=1e-4)

    # but also check that they are definitely doing something
    @test !(convF ≈ convF_5sigma)
    @test !(convF_mat ≈ convF_mat5)
end

@testset "rotation" begin
    # this implementation is less accurate than the one in Korg, but it produces correct-ish results
    function naive_apply_rotation(flux, wls::R, vsini, ε=0.6) where R <: AbstractRange
        vsini *= 1e5 # km/s to cm/s
        newFtype = promote_type(eltype(flux), eltype(wls), typeof(vsini), typeof(ε))
        newF = zeros(newFtype, length(flux))
        
        c1 = 2(1-ε)
        c2 = π * ε / 2
        
        # step(wls) makes things normalized on the grid, and the factor of v_L in Gray becomes Δλrot 
        # (because we are working in wavelenths) and moves inside the loop
        denominator = π * (1-ε/3) / step(wls) 
        
        for i in 1:length(flux)
            Δλrot = wls[i] * vsini / Korg.c_cgs
            nwls = Int(floor(Δλrot / step(wls)))
            window = max(1, i-nwls) : min(length(flux), i+nwls)
                    
            x = (wls[i] .- wls[window]) ./ Δλrot
            one_less_x2 = @. 1 - x^2
            
            @. newF[window] .+= flux[i] * (c1*sqrt(one_less_x2) + c2*one_less_x2) / (denominator * Δλrot)
        end
        newF
    end

    wls = 4090:0.01:5010
    flux = zeros(length(wls))
    flux[990:1010] .= 1

    @testset for vsini in [0.0, 1e-10, 1.0, 5.0, 10.0, 20.0], ε in [0.1, 0.6, 0.9]
        # also test handling of multiple wl ranges
        @testset for wls in [wls, [4090:0.01:5007, 5007.01:0.01:5010]]
            rflux = Korg.apply_rotation(flux, wls, vsini, ε)
            rflux2 = Korg.apply_rotation(flux, wls * 1e-8, vsini, ε)


            # rotational kernel is normalized
            @test sum(flux) ≈ sum(rflux) rtol=1e-2
            @test sum(flux) ≈ sum(rflux2) rtol=1e-2

            @test rflux == rflux2 # wl units shouldn't matter

            if vsini > 1.0 && wls isa AbstractRange
                @test assert_allclose_grid(rflux, naive_apply_rotation(flux, wls, vsini, ε), 
                                           [("λ", wls*1e8, "Å")]; atol=1e-2, print_rachet_info=false)
            end
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

@testset "linelists" begin
    atm = read_model_atmosphere("data/sun.mod")

    # iterate over fns, not lists, because they make the output of the test suite way too long
    @testset for linelist_fn in [Korg.get_VALD_solar_linelist, 
                                 Korg.get_APOGEE_DR17_linelist,
                                 Korg.get_GALAH_DR3_linelist]
        linelist = linelist_fn()
        @test issorted(linelist, by=l->l.wl)

        # make sure things run (types have caused problems in the past)
        synthesize(atm, linelist, format_A_X(), 5000, 5000)
    end
end

end #top-level testset
