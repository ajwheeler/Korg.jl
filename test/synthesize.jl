@testset "synthesize" begin
    # use this for everything
    atm_small = let atm = read_model_atmosphere("data/sun.mod")
        Korg.PlanarAtmosphere(atm.layers[40:43]) # just a few layers for fast tests
    end

    @testset "vmic specification" begin
        linelist = [Korg.Line(6000e-8, -1.0, Korg.species"Fe I", 0.3)]

        sol1 = synthesize(atm_small, linelist, format_A_X(), 6000, 6000; vmic=1.0)
        sol2 = synthesize(atm_small, linelist, format_A_X(), 6000, 6000; vmic=[1.0, 1.0, 1.0, 1.0])
        @test sol1.flux == sol2.flux

        sol3 = synthesize(atm_small, linelist, format_A_X(), 6000, 6000; vmic=[1.0, 2.0, 3.0, 4.0])
        @test sol1.flux != sol3.flux
    end

    @testset "line buffer" begin
        # strong line at 5999 Å
        line1 = Korg.Line(5999e-8, 1.0, Korg.species"Na I", 0.0)
        # strong line at 5997 Å
        line2 = Korg.Line(5997e-8, 1.0, Korg.species"Na I", 0.0)

        # use a 1 Å line buffer so only line1 in included
        sol_no_lines = synthesize(atm_small, [], format_A_X(), 6000, 6000; line_buffer=1.0) #synthesize at 6000 Å only
        sol_one_lines = synthesize(atm_small, [line1], format_A_X(), 6000, 6000; line_buffer=1.0)
        sol_two_lines = synthesize(atm_small, [line1, line2], format_A_X(), 6000, 6000;
                                   line_buffer=1.0)

        @test sol_no_lines.flux != sol_one_lines.flux
        @test sol_two_lines.flux == sol_one_lines.flux
    end

    @testset "wavelength handling" begin
        wls = Korg.Wavelengths(5000, 5001)

        msg = "The air_wavelengths keyword argument is deprecated"
        @test_warn msg synthesize(atm_small, [], format_A_X(), 5000, 5001; air_wavelengths=true)

        sol = synthesize(atm_small, [], format_A_X(), wls)
        @testset for wl_params in [
            [5000:0.01:5001],
            [[5000:0.01:5001]],
            [collect(5000:0.01:5001)]
        ]
            sol2 = synthesize(atm_small, [], format_A_X(), wl_params...)
            @test sol.flux ≈ sol2.flux
        end
    end

    @testset "MHD for H lines in APOGEE warning" begin
        msg = "if you are synthesizing at wavelengths longer than 15000 Å (e.g. for APOGEE), setting use_MHD_for_hydrogen_lines=false is recommended for the most accurate synthetic spectra. This behavior may become the default in Korg 1.0."

        @test_warn msg synthesize(atm_small, [], format_A_X(), 15000, 15001)
        @test_nowarn synthesize(atm_small, [], format_A_X(), 15000, 15001;
                                use_MHD_for_hydrogen_lines=false)
        @test_nowarn synthesize(atm_small, [], format_A_X(), 5000, 5001)
    end

    @testset "precomputed chemical equilibrium" begin
        # test that the precomputed chemical equilibrium works
        A_X = format_A_X()
        sol = synthesize(atm_small, [], A_X, 5000, 5000)
        sol_eq = synthesize(atm_small, [], A_X, 5000, 5000; use_chemical_equilibrium_from=sol)
        @test sol.flux == sol_eq.flux
    end

    @testset "mu specification" begin
        sol = synthesize(atm_small, [], format_A_X(), 5000, 5000; mu_values=5,
                         I_scheme="linear")
        @test length(sol.mu_grid) == 5

        sol = synthesize(atm_small, [], format_A_X(), 5000, 5000; mu_values=0:0.5:1.0,
                         I_scheme="linear")
        @test length(sol.mu_grid) == 3

        sol = synthesize(atm_small, [], format_A_X(), 5000, 5000; mu_values=0:0.5:1.0,
                         I_scheme="linear_flux_only")
        @test sol.mu_grid == [(1, 1)]
    end

    @testset "linelist checking" begin
        msg = "The provided linelist was not empty"
        linelist = [Korg.Line(5000e-8, 1.0, Korg.species"Na I", 0.0)]
        @test_warn msg synthesize(atm_small, linelist, format_A_X(), 6000, 6000)
    end

    @testset "linelist filtering" begin
        function naive_filter(linelist, wl_ranges, line_buffer)
            filter(linelist) do line
                map(wl_ranges) do r
                    r[1] - line_buffer <= line.wl <= r[end] + line_buffer
                end |> any
            end
        end

        b = 2.0 * 1e-8 # line_buffer
        linelist = Korg.get_APOGEE_DR17_linelist(; include_water=false)
        @testset "linelist filtering" for wls in [
            6000:7000,
            15100:15200,
            [6000:7000, 15100:15200, 16900:17100],
            [15100:15101, 15102:15103]
        ]
            wls = Korg.Wavelengths(wls)
            @test length(Korg.filter_linelist(linelist, wls, b)) ==
                  length(naive_filter(linelist, wls.wl_ranges, b))
        end
    end

    @testset "α(5000 Å) linelist" begin
        # test automatic construction of a linelist at 5000 Å for the calculation of α(5000 Å),
        # which is used by the default RT scheme

        # synthesis linelist
        ll = filter(Korg.get_VALD_solar_linelist()) do line
            4980 < line.wl * 1e8 < 5100
        end

        # if there's full coverage, don't insert anything
        @test issubset(Korg.get_alpha_5000_linelist(ll), ll)

        #if there's no coverage, use the fallback linelist
        @test Korg.get_alpha_5000_linelist([]) == Korg._alpha_5000_default_linelist

        # if there's partial coverage, insert the fallback linelist where needed
        small_ll = filter(ll) do line
            line.wl * 1e8 > 5000
        end
        ll5 = Korg.get_alpha_5000_linelist(small_ll)
        @test issorted([line.wl for line in ll5])
        # test that it transitions between the two linelists correctly
        i = findfirst(ll5 .== small_ll[1])
        @test issubset(ll5[1:i-1], Korg._alpha_5000_default_linelist)
        @test issubset(ll5[i:end], small_ll)

        small_ll = filter(ll) do line
            line.wl * 1e8 < 4995
        end
        ll5 = Korg.get_alpha_5000_linelist(small_ll)
        @test issorted([line.wl for line in ll5])
        # test that it transitions between the two linelists correctly'
        i = findfirst(ll5 .== small_ll[end])
        @test issubset(ll5[1:i], small_ll)
        @test issubset(ll5[i+1:end], Korg._alpha_5000_default_linelist)
    end

    @testset "wavelength bounds" begin
        # Test wavelength below minimum
        @test_throws ArgumentError synthesize(atm_small, [], format_A_X(), 1200, 1400)
    end
end
