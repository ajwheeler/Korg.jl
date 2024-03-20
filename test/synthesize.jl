@testset "synthesize" begin
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

    @testset "precomputed chemical equilibrium" begin
        # test that the precomputed chemical equilibrium works
        atm = read_model_atmosphere("data/sun.mod")
        A_X = format_A_X()
        sol = synthesize(atm, [], A_X, 5000, 5000)
        sol_eq = synthesize(atm, [], A_X, 5000, 5000; use_chemical_equilibrium_from=sol)
        @test sol.flux == sol_eq.flux
    end

    @testset "synthesize wavelength handling" begin
        # these are essentially tests of Korg.construct_wavelength_ranges, but we test "through"
        # synthesize because it's crucial that synthesize works

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
        sol1 = synthesize(atm, [], A_X, r1; hydrogen_lines=true)
        sol2 = synthesize(atm, [], A_X, [r2]; hydrogen_lines=true)
        sol3 = synthesize(atm, [], A_X, [r1, r2]; hydrogen_lines=true)

        @test sol1.wavelengths == sol3.wavelengths[sol3.subspectra[1]]
        @test sol2.wavelengths == sol3.wavelengths[sol3.subspectra[2]]
        @test sol1.flux == sol3.flux[sol3.subspectra[1]]
        @test sol2.flux == sol3.flux[sol3.subspectra[2]]
    end

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

    @testset "linelist checking" begin
        msg = "The provided linelist was not empty"
        atm = interpolate_marcs(5000.0, 4.4)
        linelist = [Korg.Line(5000e-8, 1.0, Korg.species"Na I", 0.0)]
        @test_warn msg synthesize(atm, linelist, format_A_X(), 6000, 6000)
    end
end