using Random

@testset "Fit" begin
    @testset "parameter scaling" begin
        params = Dict("Teff"=>3200.0, "logg"=>4.5, "m_H"=>-2.0, "vmic"=>3.2, "vsini"=>10.0, "O"=>-1.0)
        sparams = Korg.Fit.scale(params)
        uparams = Korg.Fit.unscale(sparams)
        @test all(isapprox.(values(uparams), values(params); rtol=1e-3))
    end

    @testset "fit param validation" begin
        @test_throws ArgumentError Korg.Fit.validate_params((; Teff = 3200), (;))
        @test_throws ArgumentError Korg.Fit.validate_params((; logg = 3200), (;))
        @test_throws ArgumentError Korg.Fit.validate_params((Teff=4500, logg = 3200, m_H = 0.1), (; m_H=0.1))

        # alpha may be specified, but it has no default, making it a special case
        p0, _ = Korg.Fit.validate_params((Teff=4500, logg=4.5, alpha_H=0.2), (;))
        @test p0["alpha_H"] == 0.2
        _, fixed = Korg.Fit.validate_params((Teff=4500, logg=4.5), (; alpha_H=0.2))
        @test fixed["alpha_H"] == 0.2
        p0, fixed = Korg.Fit.validate_params((Teff=4500, logg=4.5), (;))
        @test !("alpha_H" in keys(p0)) && !("alpha_H" in keys(fixed))

        for initial_guess in [(Teff=4500, logg=4.5), Dict("Teff"=> 4500, "logg"=>4.5)]
            for fixed_params in [(;), Dict()]
                _, fixed_params = Korg.Fit.validate_params(initial_guess, fixed_params)
                @test fixed_params["m_H"] == 0
                @test fixed_params["vmic"] == 1
            end
        end
    end

    @testset "merge bounds, masks etc" begin
        mbounds = Korg.Fit.merge_bounds([(1, 3), (2, 4), (5,6)], 0)
        @test mbounds[1] == [(1, 4), (5,6)]
        @test mbounds[2] == [[1, 2], [3]]

        mbounds = Korg.Fit.merge_bounds([(1, 3), (2, 6), (5,7)], 0)
        @test mbounds[1] == [(1, 7)]
        @test mbounds[2] == [[1,2,3]]

        mbounds = Korg.Fit.merge_bounds([(2, 6), (5,7), (1, 3)], 0)
        @test mbounds[1] == [(1, 7)]
        @test mbounds[2] == [[3, 1, 2]]

        obs_wls = 5000:1.0:5010
        synth_wls = 5000 : 0.01 : 5012
        windows = [(5001.0, 5002.0), (5003.0, 5004.0), (5007.0, 5008.0)]

        windows, _ = Korg.Fit.merge_bounds(windows, 1.0)
        obs_wl_mask, synth_wl_mask, multi_synth_wls = 
            Korg.Fit.calculate_multilocal_masks_and_ranges(windows, obs_wls, synth_wls, 1.0)

        @test issorted(multi_synth_wls, by=first)
        @test multi_synth_wls == [5000.0:0.01:5005.0, 5006.0:0.01:5009.0]
        @test obs_wl_mask == Bool[0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]
        @test synth_wls[synth_wl_mask] == [(5000.0:0.01:5005.0)... ;  (5006.0:0.01:5009.0)...]
    end

    @testset "fit_spectrum" begin
        # fit params
        Teff = 6402.0
        m_H = -1.02

        # fixed_params
        logg = 4.52
        vmic = 0.83
        cntm_offset = 0.04

        synth_wls = 5000:0.01:5010
        obs_wls = 5003 : 0.03 : 5008

        linelist = Korg.get_VALD_solar_linelist()

        LSF = Korg.compute_LSF_matrix(synth_wls, obs_wls, 50_000; verbose=false)

        # generate a spectrum and 
        atm = interpolate_marcs(Teff, logg, m_H)
        sol = synthesize(atm, linelist, format_A_X(m_H), [synth_wls]; vmic=vmic)
        spectrum = LSF * (sol.flux ./ (sol.cntm .* (1 - cntm_offset)))
        err = 0.01 * ones(length(spectrum)) # don't actually apply error to keep tests deterministic

        # now fit it
        p0 = (Teff=5350.0, m_H=0.0)
        fixed = (logg=logg, vmic=vmic, cntm_offset=cntm_offset)
        result = Korg.Fit.fit_spectrum(obs_wls, spectrum, err, linelist, p0, fixed; 
                                       synthesis_wls=synth_wls, LSF_matrix=LSF)
        
        params, Σ = result.covariance
        Teff_index = findfirst(params .== "Teff")
        Teff_sigma = sqrt(Σ[Teff_index, Teff_index])
        m_H_index = findfirst(params .== "m_H")
        m_H_sigma = sqrt(Σ[m_H_index, m_H_index])
        
        # check that inferred parameters are within 2 sigma of the true values
        @test result.best_fit_params["Teff"] ≈ Teff atol=1Teff_sigma
        @test result.best_fit_params["m_H"] ≈ m_H atol=1m_H_sigma

        # check that best-fit flux is within 1% of the true flux at all pixels
        @test assert_allclose(spectrum, result.best_fit_flux, rtol=0.01)

        # test argument checks
        @test_throws "must all have the same length" Korg.Fit.fit_spectrum(obs_wls[1:end-1], spectrum, err, linelist, p0, fixed; 
                                                                           synthesis_wls=synth_wls, LSF_matrix=LSF)
        @test_throws "the first dimension of LSF_matrix" Korg.Fit.fit_spectrum(obs_wls, spectrum, err, linelist, p0, fixed; 
                                                                           synthesis_wls=synth_wls, LSF_matrix=LSF[1:end-1, :])
        @test_throws "the second dimension of LSF_matrix" Korg.Fit.fit_spectrum(obs_wls, spectrum, err, linelist, p0, fixed; 
                                                                           synthesis_wls=synth_wls, LSF_matrix=LSF[:, 1:end-1])
    end

    @testset "ews_to_abundances" begin
        @testset "require sorted linelists" begin
            sun_Teff, sun_logg, sun_Fe_H, sun_vmic = (5777, 4.44, 0.0, 1.0)
            sun_A_X = Korg.format_A_X(sun_Fe_H)
            sun_atm = Korg.read_model_atmosphere("data/sun.mod")

            linelist = [
                Korg.Line(5054.642 * 1e-8, -1.92100, Korg.Species("26.0"), 3.64, 4.68e-32),
                Korg.Line(5044.211 * 1e-8, -2.05800, Korg.Species("26.0"), 2.8512, 2.71e-31),
            ]
            sun_ews = [74.3, 40.5]
            @test_throws ArgumentError Korg.Fit.ews_to_abundances(sun_atm, linelist, sun_A_X, sun_ews, vmic=sun_vmic)
        end

        @testset "length of linelist and ews should match" begin
            sun_Teff, sun_logg, sun_Fe_H, sun_vmic = (5777, 4.44, 0.0, 1.0)
            sun_A_X = Korg.format_A_X(sun_Fe_H)
            sun_atm = Korg.read_model_atmosphere("data/sun.mod")
            linelist = [Korg.Line(5044.211 * 1e-8, -2.05800, Korg.Species("26.0"), 2.8512, 2.71e-31)]
            @test_throws ArgumentError Korg.Fit.ews_to_abundances(sun_atm, linelist, sun_A_X, [], vmic=sun_vmic)
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
            @test_throws ArgumentError Korg.Fit.ews_to_abundances(sun_atm, linelist, sun_A_X, sun_ews, vmic=sun_vmic)
        end    

        @testset "Melendez et al. (2014) sanity check" begin

            sun_ews = [74.3, 40.5, 96.1, 19.1, 80.7]
            sco_ews = [74.8, 40.9, 97.5, 18.9, 84.0]
            linelist = [
                Korg.Line(5044.211 * 1e-8, -2.05800, Korg.Species("26.0"), 2.8512, 2.71e-31),
                Korg.Line(5054.642 * 1e-8, -1.92100, Korg.Species("26.0"), 3.64, 4.68e-32),
                Korg.Line(5127.359 * 1e-8, -3.30700, Korg.Species("26.0"), 0.915, 1.84e-32),
                # don't convert these ones to cm.  It should still work.
                Korg.Line(5127.679, -6.12500, Korg.Species("26.0"), 0.052, 1.2e-32),
                Korg.Line(5197.577, -2.22000, Korg.Species("26.1"), 3.2306, 8.69e-33),
            ]
            sun_Teff, sun_logg, sun_Fe_H, sun_vmic = (5777, 4.44, 0.0, 1.0)
            sun_A_X = Korg.format_A_X(sun_Fe_H)
            sun_atm = Korg.read_model_atmosphere("data/sun.mod")

            sco_teff, sco_logg, sco_fe_h, sco_vmic = (5823, 4.45, 0.054, sun_vmic + 0.02)
            sco_A_X = Korg.format_A_X(sco_fe_h)
            # Note: NOT true, but just for testing the whole performance
            sco_atm = sun_atm
            
            sun_abundances = Korg.Fit.ews_to_abundances(sun_atm, linelist, sun_A_X, sun_ews, vmic=sun_vmic)
            sco_abundances = Korg.Fit.ews_to_abundances(sco_atm, linelist, sco_A_X, sco_ews, vmic=sco_vmic)        
            diff_abundances = sco_abundances .- sun_abundances

            mean_diff_abundances = sum(diff_abundances) / length(diff_abundances)
            @test abs(mean_diff_abundances) < 0.05
            # TODO: test for stddev?        
        end
    end

    @testset "stellar parameters via EWs" begin
        # used in both tests below
        good_linelist = [Korg.Line(5e-5, -2.05800, Korg.species"Fe I",  4.1),
                 Korg.Line(6e-5, -1.92100, Korg.species"Fe I",  4.2),
                 Korg.Line(7e-5, -1.92100, Korg.species"Fe I",  4.3),
                 Korg.Line(8e-5, -1.92100, Korg.species"Fe II", 4.4)]

        @testset "validate arguments" begin
            # vmic can't be specified
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(good_linelist, ones(length(good_linelist)); vmic=1.0)

            # number of EWs, EW uncertainties, and lines should match
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(good_linelist, [1.0])
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(good_linelist, ones(length(good_linelist)), [1.0])

            # different elements
            linelist = [Korg.Line(5e-5, -2.05800, Korg.species"Fe I", 0, 0),
                        Korg.Line(6e-5, -1.92100, Korg.species"Fe I", 0, 0),
                        Korg.Line(7e-5, -1.92100, Korg.species"Mn I", 0, 0),
                        Korg.Line(8e-5, -1.92100, Korg.species"Fe II", 0, 0)]
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(linelist, ones(length(linelist)))

            # can't use a molecule
            linelist = [Korg.Line(5e-5, -2.05800, Korg.species"CO I", 0, 0),
                        Korg.Line(6e-5, -1.92100, Korg.species"CO I", 0, 0),
                        Korg.Line(7e-5, -1.92100, Korg.species"CO I", 0, 0),
                        Korg.Line(8e-5, -1.92100, Korg.species"CO II", 0, 0)]
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(linelist, ones(length(linelist)))

            # no ions
            linelist = [Korg.Line(5e-5, -2.05800, Korg.species"Fe I", 0, 0),
                        Korg.Line(6e-5, -1.92100, Korg.species"Fe I", 0, 0),
                        Korg.Line(6.5e-5, -1.92100, Korg.species"Fe I", 0, 0),
                        Korg.Line(7e-5, -1.92100, Korg.species"Fe I", 0, 0)]
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(linelist, ones(length(linelist)))

            # not enough lines
            linelist = [Korg.Line(5e-5, -2.05800, Korg.species"Fe I", 0, 0),
                        Korg.Line(6e-5, -1.92100, Korg.species"Fe I", 0, 0),
                        Korg.Line(7e-5, -1.92100, Korg.species"Fe II", 0, 0)]
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(linelist, ones(length(linelist)))

            # parameter ranges must have positive measure
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(good_linelist, ones(length(good_linelist)); 
                        parameter_ranges=[(5777, 5777), (3.0, 4.0), (1.0, 2.0), (-0.5, 0.0)])
            
            # vmic0 can't be 0
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(good_linelist, ones(length(good_linelist)); vmic0=0)
        end

        @testset "basic fit" begin
            # 2 Å wide window around each line
            synth_wls = map(good_linelist) do line
                wl = line.wl * 1e8
                wl - 1.0 : 0.01 : wl + 1.0
            end
            sol = synthesize(interpolate_marcs(5777.0, 4.44), good_linelist, format_A_X(),
                             synth_wls; hydrogen_lines=false)

            # EWs you get for the fake linelist with solar params
            # the real implementation uses the trapezoid rule, but this is close enough
            EWs = [sum((1 .- sol.flux ./ sol.cntm)[r]) for r in sol.subspectra] * 10 #0.01 Å -> mÅ
            EW_err = ones(length(EWs)) * 0.5
            best_fit_params, stat_err, sys_err = Korg.Fit.ews_to_stellar_parameters(good_linelist, EWs, EW_err)

            @test sys_err == [0.0, 0.0, 0.0, 0.0]
            for (i, p) in enumerate([5777.0, 4.44, 1.0, 0.0])
                @test best_fit_params[i] ≈ p atol=stat_err[i]
            end

            # check that fixing parameters works
            for (i, param) in enumerate([:Teff0, :logg0, :vmic0, :m_H0])
                fixed_params = zeros(Bool, 4)
                fixed_params[i] = true
                kwargs = Dict(param => best_fit_params[i])
                # start teff and logg close to right answer to make it faster.
                bestfit_fixed, _, _ = Korg.Fit.ews_to_stellar_parameters(
                    good_linelist, EWs, EW_err; Teff0=5740.0, logg0=4.4, fix_params=fixed_params, kwargs...)

                for i in 1:4
                    @test bestfit_fixed[i] ≈ best_fit_params[i] atol=stat_err[i]/10
                end
            end
        end
    end
end