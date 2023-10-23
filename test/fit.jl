@testset "Fit" begin
    @testset "parameter scaling" begin
        params = (Teff = 3200, logg=4.5, m_H=-2, vmic=3.2, vsini=10, O=-1)
        sparams = Korg.Fit.scale(params)
        uparams = Korg.Fit.unscale(sparams)
        @test all(isapprox.(values(uparams), values(params); rtol=1e-3))
    end

    @testset "fit param validation" begin
        @test_throws ArgumentError Korg.Fit.validate_params((; Teff = 3200), (;))
        @test_throws ArgumentError Korg.Fit.validate_params((; logg = 3200), (;))
        @test_throws ArgumentError Korg.Fit.validate_params((Teff=4500, logg = 3200, m_H = 0.1), (; m_H=0.1))

        _, fixed_params = Korg.Fit.validate_params((Teff=4500, logg=4.5), (;))
        @test fixed_params.m_H == 0
        @test fixed_params.vmic == 1
    end

    @testset "merge bounds, masks etc" begin
        @test Korg.Fit.merge_bounds([(1, 3), (2, 4), (5,6)], 0) == [(1, 4), (5,6)]
        @test Korg.Fit.merge_bounds([(1, 3), (2, 6), (5,7)], 1.0) == [(1, 7)]

        obs_wls = 5000:1.0:5010
        synth_wls = 5000 : 0.01 : 5012
        windows = [(5001.0, 5002.0), (5003.0, 5004.0), (5007.0, 5008.0)]

        windows = Korg.Fit.merge_bounds(windows, 1.0)
        obs_wl_mask, synth_wl_mask, multi_synth_wls = 
            Korg.Fit.calculate_multilocal_masks_and_ranges(windows, obs_wls, synth_wls, 1.0)

        @test issorted(multi_synth_wls, by=first)
        @test multi_synth_wls == [5000.0:0.01:5005.0, 5006.0:0.01:5009.0]
        @test obs_wl_mask == Bool[0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]
        @test synth_wls[synth_wl_mask] == [(5000.0:0.01:5005.0)... ;  (5006.0:0.01:5009.0)...]
    end

    @testset "don't allow hydrogen lines in ew_to_abundances" begin
        sun_Teff, sun_logg, sun_Fe_H, sun_vmic = (5777, 4.44, 0.0, 1.0)
        sun_A_X = Korg.format_A_X(sun_Fe_H)
        sun_atm = Korg.read_model_atmosphere("data/sun.mod")

        linelist = [Korg.Line(5044.211 * 1e-8, -2.05800, Korg.Species("26.0"), 2.8512, 2.71e-31)]
        sun_ews = [74.3]
        @test_throws ArgumentError Korg.ews_to_abundances(sun_atm, linelist, sun_A_X, sun_ews, vmic=sun_vmic, hydrogen_lines=true)        
    end
    @testset "require sorted linelists" begin
        sun_Teff, sun_logg, sun_Fe_H, sun_vmic = (5777, 4.44, 0.0, 1.0)
        sun_A_X = Korg.format_A_X(sun_Fe_H)
        sun_atm = Korg.read_model_atmosphere("data/sun.mod")

        linelist = [
            Korg.Line(5054.642 * 1e-8, -1.92100, Korg.Species("26.0"), 3.64, 4.68e-32),
            Korg.Line(5044.211 * 1e-8, -2.05800, Korg.Species("26.0"), 2.8512, 2.71e-31),
        ]
        sun_ews = [74.3, 40.5]
        @test_throws ArgumentError Korg.ews_to_abundances(sun_atm, linelist, sun_A_X, sun_ews, vmic=sun_vmic)
    end

    @testset "length of linelist and ews should match" begin
        sun_Teff, sun_logg, sun_Fe_H, sun_vmic = (5777, 4.44, 0.0, 1.0)
        sun_A_X = Korg.format_A_X(sun_Fe_H)
        sun_atm = Korg.read_model_atmosphere("data/sun.mod")
        linelist = [Korg.Line(5044.211 * 1e-8, -2.05800, Korg.Species("26.0"), 2.8512, 2.71e-31)]
        @test_throws ArgumentError Korg.ews_to_abundances(sun_atm, linelist, sun_A_X, [], vmic=sun_vmic)
    end

    @testset "disallow molecules" begin
        sun_Teff, sun_logg, sun_Fe_H, sun_vmic = (5777, 4.44, 0.0, 1.0)
        sun_A_X = Korg.format_A_X(sun_Fe_H)
        sun_atm = Korg.read_model_atmosphere("data/sun.mod")
        linelist = [
            Korg.Line(5044.211 * 1e-8, -2.05800, Korg.Species("26.0"), 2.8512, 2.71e-31),
            Korg.Line(5044.211 * 1e-8, -2.05800, Korg.Species("106.0"), 2.8512, 2.71e-31)

        ]
        sun_ews = [74.3, 40.5]
        @test_throws ArgumentError Korg.ews_to_abundances(sun_atm, linelist, sun_A_X, sun_ews, vmic=sun_vmic)
    end    

    @testset "Melendez et al. (2014) sanity check" begin

        sun_ews = [74.3, 40.5, 96.1, 19.1, 80.7]
        sco_ews = [74.8, 40.9, 97.5, 18.9, 84.0]
        linelist = [
            Korg.Line(5044.211 * 1e-8, -2.05800, Korg.Species("26.0"), 2.8512, 2.71e-31),
            Korg.Line(5054.642 * 1e-8, -1.92100, Korg.Species("26.0"), 3.64, 4.68e-32),
            Korg.Line(5127.359 * 1e-8, -3.30700, Korg.Species("26.0"), 0.915, 1.84e-32),
            Korg.Line(5127.679 * 1e-8, -6.12500, Korg.Species("26.0"), 0.052, 1.2e-32),
            Korg.Line(5197.577 * 1e-8, -2.22000, Korg.Species("26.1"), 3.2306, 8.69e-33),
        ]
        sun_Teff, sun_logg, sun_Fe_H, sun_vmic = (5777, 4.44, 0.0, 1.0)
        sun_A_X = Korg.format_A_X(sun_Fe_H)
        sun_atm = Korg.read_model_atmosphere("data/sun.mod")

        sco_teff, sco_logg, sco_fe_h, sco_vmic = (5823, 4.45, 0.054, sun_vmic + 0.02)
        sco_A_X = Korg.format_A_X(sco_fe_h)
        # Note: NOT true, but just for testing the whole performance
        sco_atm = sun_atm
        
        sun_abundances = ews_to_abundances(sun_atm, linelist, sun_A_X, sun_ews, vmic=sun_vmic)
        sco_abundances = ews_to_abundances(sco_atm, linelist, sco_A_X, sco_ews, vmic=sco_vmic)        
        diff_abundances = sco_abundances .- sun_abundances

        mean_diff_abundances = sum(diff_abundances) / length(diff_abundances)
        @test abs(mean_diff_abundances) < 0.05
        # TODO: test for stddev?        
    end

    @testset "test neighbourhood grouping" begin
        linelist = [
            Korg.Line(5044.211 * 1e-8, -2.05800, Korg.Species("26.0"), 2.8512, 2.71e-31),
            Korg.Line(5054.642 * 1e-8, -1.92100, Korg.Species("26.0"), 3.64, 4.68e-32),
            Korg.Line(5127.359 * 1e-8, -3.30700, Korg.Species("26.0"), 0.915, 1.84e-32),
            Korg.Line(5127.679 * 1e-8, -6.12500, Korg.Species("26.0"), 0.052, 1.2e-32),
            Korg.Line(5197.577 * 1e-8, -2.22000, Korg.Species("26.1"), 3.2306, 8.69e-33),
        ]
        
        @test length(linelist_neighbourhood_indices(linelist, 2)) == 2
        @test length(linelist_neighbourhood_indices(linelist, 10)) == 4 # (5044, 5054, 5127.3, (5127.6, 5197.5))
        @test length(linelist_neighbourhood_indices(linelist[1:3], 2)) == 1
    end

end