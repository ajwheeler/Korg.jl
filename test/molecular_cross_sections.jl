@testset "precomputed molecular opacities" begin
    # use 10 Å of the APOGEE linelist for testing
    linelist = Korg.get_APOGEE_DR17_linelist()
    filter!(linelist) do line
        15020 < line.wl * 1e8 < 15030
    end

    # split into water and non-water lines
    water_mask = [l.species == Korg.species"H2O" for l in linelist]
    linelist_less_water = linelist[.!water_mask]
    water_lines = linelist[water_mask]

    wls = [15020:0.01:15025, 15030:0.01:15035]

    # a couple vmic kernels broader
    table_wls = [15019:0.01:15026, 15029:0.01:15036]

    table = Korg.MolecularCrossSection(water_lines, table_wls)
    @test !Korg.is_vmic_gridded(table) # not legacy version

    # this test fails if I truncate the atmosphere to 40:43, which _may_ indicate that it's worth
    # better understanding the precision of interpolation here
    atm = interpolate_marcs(4000.0, 4.5)

    filename = tempname()
    Korg.save_molecular_cross_section(filename, table)
    deserialized_table = Korg.read_molecular_cross_section(filename)
    @test !Korg.is_vmic_gridded(deserialized_table) # not legacy version
    @test table.itp.itp.coefs == deserialized_table.itp.itp.coefs
    @test all(table.wls .≈ deserialized_table.wls)
    @test table.species == deserialized_table.species

    depth_dependent_vmics = fill(0.5, length(atm.layers))
    depth_dependent_vmics[20:end] .= 0.0
    @testset for vmic in [0.5, 1.2, depth_dependent_vmics]
        sol_without = synthesize(atm, linelist, format_A_X(), wls; vmic=vmic)
        sol_with = synthesize(atm, linelist_less_water, format_A_X(), wls;
                              molecular_cross_sections=[deserialized_table], vmic=vmic)

        @test assert_allclose(sol_without.flux, sol_with.flux; rtol=1e-3)
    end

    @testset "vmic differentiability" begin
        Ts = [3500.0, 4500.0]

        number_densities = Dict(Korg.species"H2O" => [1e10, 1e10])
        function convolved_alpha(vmic, wavelengths)
            α = zeros(promote_type(typeof(vmic), Float64), length(Ts), length(wavelengths))
            Korg.interpolate_molecular_cross_sections!(α, [table], wavelengths, Ts, vmic,
                                                       number_densities)
            sum(α) # arbitrary, just to give the function a scalar output
        end

        # (single wavelength is like the reference wavelength in anchored radiative transfer)
        for wavelengths in Korg.Wavelengths.([wls, [15022.0]])
            function f(vmic)
                α = zeros(promote_type(typeof(vmic), Float64), length(Ts), length(wavelengths))
                Korg.interpolate_molecular_cross_sections!(α, [table], wavelengths, Ts, vmic,
                                                           number_densities)
                sum(α) # arbitrary, just to give the function a scalar output
            end

            f = vmic -> convolved_alpha(vmic, wavelengths)
            @test ForwardDiff.derivative(f, 1.0)≈
            FiniteDiff.finite_difference_derivative(f, 1.0) rtol=1e-6

            # the derivative with respect to vmic vanishes at vmic = 0
            @test ForwardDiff.derivative(f, 0.0) == 0
            # for kernels narrower than the table spacing (vmic below ~0.04 km/s here), the
            # cross-section is interpolated without broadening and the derivative is taken as 0
            @test f(0.01) == f(0.0)
            @test ForwardDiff.derivative(f, 0.01) == 0
        end
    end

    @testset "tables with too-coarse wl sampling error" begin
        coarse_table = Korg.MolecularCrossSection(water_lines, 15019:0.02:15026)
        Ts = [3500.0, 4500.0]
        number_densities = Dict(Korg.species"H2O" => [1e10, 1e10])
        synth_wls = Korg.Wavelengths(15020:0.01:15025)
        α = zeros(length(Ts), length(synth_wls))
        @test_throws ArgumentError Korg.interpolate_molecular_cross_sections!(α, [coarse_table],
                                                                              synth_wls, Ts, 1.0,
                                                                              number_densities)

        # single-wavelength requests have no spacing, so any table is acceptable
        α_ref = zeros(length(Ts), 1)
        Korg.interpolate_molecular_cross_sections!(α_ref, [coarse_table],
                                                   Korg.Wavelengths([15022.0]), Ts, 1.0,
                                                   number_densities)
        @test all(α_ref .> 0)
    end

    @testset "legacy vmic-gridded tables" begin
        legacy_table = Korg.MolecularCrossSection(water_lines, table_wls;
                                                  vmic_vals=[0.0, 0.5, 1.0])
        @test Korg.is_vmic_gridded(legacy_table)

        tmp = tempname()
        Korg.save_molecular_cross_section(tmp, legacy_table)
        deserialized_legacy = Korg.read_molecular_cross_section(tmp)

        @test Korg.is_vmic_gridded(deserialized_legacy)
        @test legacy_table.itp.itp.coefs == deserialized_legacy.itp.itp.coefs

        sol_without = synthesize(atm, linelist, format_A_X(), wls; vmic=0.5)
        sol_legacy = synthesize(atm, linelist_less_water, format_A_X(), wls;
                                molecular_cross_sections=[deserialized_legacy], vmic=0.5)
        @test assert_allclose(sol_without.flux, sol_legacy.flux; rtol=1e-3)
    end
end
