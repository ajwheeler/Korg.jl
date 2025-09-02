using Random

# print the timing info because this is kinda slow
@testset "Fit" verbose=true showtiming=true begin
    @testset "fit_spectrum" begin
        @testset "parameter scaling" begin
            params = Dict("Teff" => 3200.0, "logg" => 4.5, "M_H" => -2.0, "vmic" => 3.2,
                          "vsini" => 10.0, "O" => -1.0)
            sparams = Korg.Fit.scale(params)
            uparams = Korg.Fit.unscale(sparams)
            @test all(isapprox.(values(uparams), values(params); rtol=1e-3))
        end

        @testset "fit param validation" begin
            @test_throws ArgumentError Korg.Fit.validate_params((; Teff=3200), (;))
            @test_throws ArgumentError Korg.Fit.validate_params((; logg=3200), (;))
            @test_throws ArgumentError Korg.Fit.validate_params((Teff=4500, logg=3200, M_H=0.1),
                                                                (; M_H=0.1))

            @test_throws "Use M_H instead" Korg.Fit.validate_params((; m_H=0.1, Teff=5000,
                                                                     logg=4.5), (;))
            @test_throws "Use M_H instead" Korg.Fit.validate_params((; Teff=5000, logg=4.5),
                                                                    (; m_H=0.1))

            # alpha may be specified, but it has no default, making it a special case
            p0, _ = Korg.Fit.validate_params((Teff=4500, logg=4.5, alpha_H=0.2), (;))
            @test p0["alpha_H"] == 0.2
            _, fixed = Korg.Fit.validate_params((Teff=4500, logg=4.5), (; alpha_H=0.2))
            @test fixed["alpha_H"] == 0.2
            p0, fixed = Korg.Fit.validate_params((Teff=4500, logg=4.5), (;))
            @test !("alpha_H" in keys(p0)) && !("alpha_H" in keys(fixed))

            for initial_guess in [(Teff=4500, logg=4.5), Dict("Teff" => 4500, "logg" => 4.5)]
                for fixed_params in [(;), Dict()]
                    _, fixed_params = Korg.Fit.validate_params(initial_guess, fixed_params)
                    @test fixed_params["M_H"] == 0
                    @test fixed_params["vmic"] == 1
                end
            end

            @test_warn "Instead of using the `\"cntm_offset\"`` and `\"cntm_slope\"` " begin
                Korg.Fit.validate_params((; Teff=4500, cntm_slope=0), (; logg=4.0))
            end
            @test_warn "Instead of using the `\"cntm_offset\"`` and `\"cntm_slope\"` " begin
                Korg.Fit.validate_params((; Teff=4500, cntm_offset=0), (; logg=4.0))
            end
        end

        # it would be good to write some tests for Korg.Fit._setup_wavelengths_and_LSF directly

        @testset "calls to fit_spectrum" begin
            perturb!(flux, _, _) = (flux .= flux .^ 1.2)
            linelist = Korg.get_VALD_solar_linelist()

            obs_wls = 4997:0.07:5010
            windows = [(5000, 5003), (5008, 5010)]

            # make these a superset of the windows
            synth_wls = Korg.Wavelengths((4999, 5011))

            R = 50_000
            LSF = Korg.compute_LSF_matrix(synth_wls, obs_wls, R)

            fixed_params = (; vmic=0.83, logg=4.52)
            # true parameters
            Teff = 5350.0
            M_H = 0.11

            # start slightly off.  It would be a more robust test to start further away, but it
            # would also be more expensive.
            p0 = (; Teff=5000.0, M_H=0.0)

            # synthesize a spectrum, to be turned into fake data
            A_X = Korg.format_A_X(M_H, M_H)
            atm = interpolate_marcs(Teff, fixed_params.logg, A_X)
            sol = synthesize(atm, linelist, A_X, synth_wls; vmic=fixed_params.vmic)

            fake_data = LSF * (sol.flux ./ sol.cntm)
            # add an "instrumental" effect to be taken out with postprocess
            fake_data .= fake_data .^ 1.2

            # apply LSF and resampling, and put a crazy continuum in each window
            # it's important to do this after perturbing the data.
            for (i, (λstart, λstop)) in enumerate(windows)
                m = λstart .<= obs_wls .<= λstop

                slope = i * 0.01
                offset = 1.1 - 5000 * slope

                fake_data[m] .*= offset .+ slope * obs_wls[m]
            end

            # bad pixels outside the windows should be OK
            fake_data[1] = NaN
            fake_data[2] = Inf

            err = 0.01 * ones(length(obs_wls)) # don't actually apply error to keep tests deterministic

            @testset "fit test" begin
                result = Korg.Fit.fit_spectrum(obs_wls, fake_data, err, linelist, p0, fixed_params;
                                               precision=1e-2, # loose tolerances for speed
                                               postprocess=perturb!,
                                               adjust_continuum=true,
                                               R=R,
                                               windows=windows,)

                params, Σ = result.covariance

                Teff_index = findfirst(params .== "Teff")
                Teff_sigma = sqrt(Σ[Teff_index, Teff_index])
                M_H_index = findfirst(params .== "M_H")
                M_H_sigma = sqrt(Σ[M_H_index, M_H_index])

                # check that inferred parameters are within 1 sigma of the true values
                @test result.best_fit_params["Teff"]≈Teff atol=Teff_sigma
                @test result.best_fit_params["M_H"]≈M_H atol=M_H_sigma

                # check that best-fit flux is within 1% of the true flux at all pixels
                @test assert_allclose(fake_data[result.obs_wl_mask], result.best_fit_flux,
                                      rtol=0.001)
            end

            # get rid of the bad pixels to make it easer to test the error messages below
            fake_data[1:3] .= 0

            @testset "argument checks" begin
                @test_throws "LSF_matrix and synthesis_wls cannot be specified if R is provided" begin
                    Korg.Fit.fit_spectrum(obs_wls, fake_data, err, linelist, p0, fixed_params;
                                          R=R, synthesis_wls=synth_wls, LSF_matrix=LSF,
                                          time_limit=1)
                end
                @test_throws "LSF_matrix and synthesis_wls can't be specified if windows is also specified" begin
                    Korg.Fit.fit_spectrum(obs_wls, fake_data, err, linelist, p0, fixed_params;
                                          synthesis_wls=synth_wls, LSF_matrix=LSF,
                                          windows=[(5003, 5008)], time_limit=1)
                end
                @test_throws "obs_wls, obs_flux, and obs_err must all have the same length" begin
                    Korg.Fit.fit_spectrum(obs_wls[1:end-1], fake_data, err, linelist, p0,
                                          fixed_params;
                                          synthesis_wls=synth_wls, LSF_matrix=LSF[1:end-1, :],
                                          time_limit=1)
                end

                @test_throws "obs_wls must be sorted in order of increasing wavelength" begin
                    Korg.Fit.fit_spectrum(reverse(obs_wls), fake_data, err, linelist, p0,
                                          fixed_params;
                                          synthesis_wls=synth_wls, windows=windows, R=Inf,
                                          time_limit=1)
                end

                @testset "bad vals" for bad_val in [NaN, Inf, -Inf]
                    @test_throws "obs_wls, obs_flux, and obs_err must not contain NaN or Inf" begin
                        bad_flux = copy(fake_data)
                        bad_flux[1] = bad_val
                        Korg.Fit.fit_spectrum(obs_wls, bad_flux, err, linelist, p0,
                                              fixed_params;
                                              synthesis_wls=synth_wls, LSF_matrix=LSF,
                                              time_limit=1)
                    end
                end

                @testset "err contains zeros" begin
                    @test_throws "obs_err must not contain zeros" begin
                        bad_err = copy(err)
                        bad_err[1] = 0
                        Korg.Fit.fit_spectrum(obs_wls, fake_data, bad_err, linelist, p0,
                                              fixed_params;
                                              synthesis_wls=synth_wls, LSF_matrix=LSF,
                                              time_limit=1)
                    end
                end

                @test_throws "the first dimension of LSF_matrix" begin
                    Korg.Fit.fit_spectrum(obs_wls, fake_data, err, linelist, p0, fixed_params;
                                          synthesis_wls=synth_wls, LSF_matrix=LSF[1:end-1, :],
                                          time_limit=1)
                end
                @test_throws "the second dimension of LSF_matrix" begin
                    Korg.Fit.fit_spectrum(obs_wls, fake_data, err, linelist, p0, fixed_params;
                                          synthesis_wls=synth_wls, LSF_matrix=LSF[:, 1:end-1],
                                          time_limit=1)
                end
            end
        end
    end

    @testset "ews_to_abundances" begin
        @testset "require sorted linelists" begin
            sun_Teff, sun_logg, sun_Fe_H, sun_vmic = (5777, 4.44, 0.0, 1.0)
            sun_A_X = Korg.format_A_X(sun_Fe_H)
            sun_atm = Korg.read_model_atmosphere("data/sun.mod")

            linelist = [
                Korg.Line(5054.642 * 1e-8, -1.92100, Korg.Species("26.0"), 3.64, 4.68e-32),
                Korg.Line(5044.211 * 1e-8, -2.05800, Korg.Species("26.0"), 2.8512, 2.71e-31)
            ]
            sun_ews = [74.3, 40.5]
            @test_throws ArgumentError Korg.Fit.ews_to_abundances(sun_atm, linelist, sun_A_X,
                                                                  sun_ews, vmic=sun_vmic)
        end

        @testset "length of linelist and ews should match" begin
            sun_Teff, sun_logg, sun_Fe_H, sun_vmic = (5777, 4.44, 0.0, 1.0)
            sun_A_X = Korg.format_A_X(sun_Fe_H)
            sun_atm = Korg.read_model_atmosphere("data/sun.mod")
            linelist = [Korg.Line(5044.211 * 1e-8, -2.05800, Korg.Species("26.0"), 2.8512,
                                  2.71e-31)]
            @test_throws ArgumentError Korg.Fit.ews_to_abundances(sun_atm, linelist, sun_A_X,
                                                                  [],
                                                                  vmic=sun_vmic)
        end

        @testset "disallow molecules" begin
            sun_Teff, sun_logg, sun_Fe_H, sun_vmic = (5777, 4.44, 0.0, 1.0)
            sun_A_X = Korg.format_A_X(sun_Fe_H)
            sun_atm = Korg.read_model_atmosphere("data/sun.mod")
            linelist = [
                Korg.Line(5044.211 * 1e-8, -2.05800, Korg.Species("26.0"), 2.8512, 2.71e-31),
                Korg.Line(5044.211 * 1e-8, -2.05800, Korg.Species("106.0"), 2.8512, 2.71e-31)]
            sun_ews = [74.3, 40.5]
            @test_throws ArgumentError Korg.Fit.ews_to_abundances(sun_atm, linelist, sun_A_X,
                                                                  sun_ews, vmic=sun_vmic)
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
                Korg.Line(5197.577, -2.22000, Korg.Species("26.1"), 3.2306, 8.69e-33)
            ]
            sun_Teff, sun_logg, sun_Fe_H, sun_vmic = (5777, 4.44, 0.0, 1.0)
            sun_A_X = Korg.format_A_X(sun_Fe_H)
            sun_atm = Korg.read_model_atmosphere("data/sun.mod")

            sco_teff, sco_logg, sco_fe_h, sco_vmic = (5823, 4.45, 0.054, sun_vmic + 0.02)
            sco_A_X = Korg.format_A_X(sco_fe_h)
            # Note: NOT true, but just for testing the whole performance
            sco_atm = sun_atm

            sun_abundances = Korg.Fit.ews_to_abundances(sun_atm, linelist, sun_A_X, sun_ews;
                                                        vmic=sun_vmic)
            sco_abundances = Korg.Fit.ews_to_abundances(sco_atm, linelist, sco_A_X, sco_ews;
                                                        vmic=sco_vmic)
            diff_abundances = sco_abundances .- sun_abundances

            mean_diff_abundances = sum(diff_abundances) / length(diff_abundances)
            @test abs(mean_diff_abundances) < 0.05
            # consider testing stddev?
        end
    end

    @testset "stellar parameters via EWs" begin
        # used in both tests below
        simple_linelist = [Korg.Line(5e-5, -2.05800, Korg.species"Fe I", 4.1),
            Korg.Line(6e-5, -1.92100, Korg.species"Fe I", 4.2),
            Korg.Line(7e-5, -1.92100, Korg.species"Fe I", 4.3),
            Korg.Line(8e-5, -1.92100, Korg.species"Fe II", 4.4)]

        @testset "validate arguments" begin
            # vmic can't be specified
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(simple_linelist,
                                                                          ones(length(simple_linelist));
                                                                          vmic=1.0)

            # number of EWs, EW uncertainties, and lines should match
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(simple_linelist, [1.0])

            # different elements
            linelist = [Korg.Line(5e-5, -2.05800, Korg.species"Fe I", 0, 0),
                Korg.Line(6e-5, -1.92100, Korg.species"Fe I", 0, 0),
                Korg.Line(7e-5, -1.92100, Korg.species"Mn I", 0, 0),
                Korg.Line(8e-5, -1.92100, Korg.species"Fe II", 0, 0)]
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(linelist,
                                                                          ones(length(linelist)))

            # can't use a molecule
            linelist = [Korg.Line(5e-5, -2.05800, Korg.species"CO I", 0, 0),
                Korg.Line(6e-5, -1.92100, Korg.species"CO I", 0, 0),
                Korg.Line(7e-5, -1.92100, Korg.species"CO I", 0, 0),
                Korg.Line(8e-5, -1.92100, Korg.species"CO II", 0, 0)]
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(linelist,
                                                                          ones(length(linelist)))

            # no ions
            linelist = [Korg.Line(5e-5, -2.05800, Korg.species"Fe I", 0, 0),
                Korg.Line(6e-5, -1.92100, Korg.species"Fe I", 0, 0),
                Korg.Line(6.5e-5, -1.92100, Korg.species"Fe I", 0, 0),
                Korg.Line(7e-5, -1.92100, Korg.species"Fe I", 0, 0)]
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(linelist,
                                                                          ones(length(linelist)))

            # not enough lines
            linelist = [Korg.Line(5e-5, -2.05800, Korg.species"Fe I", 0, 0),
                Korg.Line(6e-5, -1.92100, Korg.species"Fe I", 0, 0),
                Korg.Line(7e-5, -1.92100, Korg.species"Fe II", 0, 0)]
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(linelist,
                                                                          ones(length(linelist)))

            # parameter ranges must have positive measure
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(simple_linelist,
                                                                          ones(length(simple_linelist));
                                                                          parameter_ranges=[(5777,
                                                                                             5777),
                                                                              (3.0, 4.0),
                                                                              (1.0, 2.0),
                                                                              (-0.5, 0.0)])

            # vmic0 can't be 0
            @test_throws ArgumentError Korg.Fit.ews_to_stellar_parameters(simple_linelist,
                                                                          ones(length(simple_linelist));
                                                                          vmic0=0)
        end

        @testset "basic fit" verbose=true showtiming=true begin
            linelist = Korg.get_VALD_solar_linelist()
            i = 0
            filter!(linelist) do line
                i += 1
                p = (Korg.get_atoms(line.species) == [26]) && (line.wl * 1e8 .> 8500)
                # take only every 10th neutral line to make it faster
                p &= (line.species.charge > 0 || i % 20 == 0)
            end

            # 2 Å wide window around each line
            synth_wls = map(linelist) do line
                wl = line.wl * 1e8
                wl-1.0:0.01:wl+1.0
            end
            A_X = format_A_X()
            sol = synthesize(interpolate_marcs(5777.0, 4.44, A_X), linelist, A_X, synth_wls;
                             hydrogen_lines=false)

            param_tolerance = [10, 0.1, 0.1, 0.01]

            # EWs you get for the fake linelist with solar params
            # the real implementation uses the trapezoid rule, but this is close enough
            EWs = [sum((1 .- sol.flux./sol.cntm)[r]) for r in sol.subspectra] * 10 #0.01 Å -> mÅ
            best_fit_params, param_err = Korg.Fit.ews_to_stellar_parameters(linelist, EWs)

            # in this fake data case, the reported uncertainties are meaningless because there is no
            # line-to-line scatter, save for what arises from imperfect fit convergence.
            # so, just check that the fit is good.
            for (guess, p, tol) in zip(best_fit_params, [5777.0, 4.44, 1.0, 0.0], param_tolerance)
                @test guess≈p atol=tol
            end

            # check that fixing parameters works
            @testset "fixed parameters" for (i, param) in enumerate([:Teff0, :logg0, :vmic0, :M_H0])
                fixed_params = zeros(Bool, 4)
                fixed_params[i] = true
                kwargs = Dict(param => best_fit_params[i])
                # start teff and logg close to right answer to make it faster.
                bestfit_fixed, param_err = Korg.Fit.ews_to_stellar_parameters(linelist, EWs;
                                                                              Teff0=5740.0,
                                                                              logg0=4.4,
                                                                              fix_params=fixed_params,
                                                                              kwargs...)

                # check that the fixed parameters are close to the one from the full fit with 10x
                # smaller tolerance
                for (fixed, best, tol) in zip(bestfit_fixed, best_fit_params, param_tolerance)
                    @test fixed≈best atol=tol/5
                end
            end
        end

        @testset "unconverged line behavior" begin
            A_X = format_A_X()
            atm = interpolate_marcs(5000, 4.0)

            # make a linelist with one too-strong line
            linelist = [Korg.Line(5000e-8, -2.5, Korg.species"Fe I", 4.0),
                Korg.Line(5001e-8, -2.5, Korg.species"Fe I", 4.0),
                Korg.Line(5002e-8, -2.5, Korg.species"Fe I", 4.0),
                Korg.Line(5003e-8, -2.5, Korg.species"Fe II", 4.0)]

            measured_EWs = Korg.Fit.calculate_EWs(atm, linelist, A_X)
            As = Korg.Fit.ews_to_abundances(atm, linelist, A_X, measured_EWs)

            # do a simple sanity check that calculate_EWs
            # and ews_to_abundances do indeed invert each other
            @test all(As .== Korg.default_solar_abundances[26])

            # make sure one of the lines fails
            erroneous_EWs = copy(measured_EWs)
            erroneous_EWs[1:2] .= 0 # mess up some neutral line measurements
            m = "Less than 70% of the lines converged"
            @test_throws m Korg.Fit.ews_to_stellar_parameters(linelist, erroneous_EWs)

            erroneous_EWs = copy(measured_EWs)
            erroneous_EWs[end] = 0 # mess up the one ion line measurement
            m = "Less than 50% of the ion lines converged"
            @test_throws m Korg.Fit.ews_to_stellar_parameters(linelist, erroneous_EWs)
        end
    end
end
