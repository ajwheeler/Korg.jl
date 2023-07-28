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
end