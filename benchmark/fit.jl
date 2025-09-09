let # for namespace hygiene
    ll, ul = 5000, 5100
    synth_wls = Korg.Wavelengths((ll, ul))
    obs_wls = (ll-10):0.07:(ul+10)

    linelist = filter(Korg.get_VALD_solar_linelist()) do line
        synth_wls[1] < line.wl < synth_wls[2] # each of these is in cm
    end

    R = 50_000

    # true parameters
    Teff = 5350.0
    M_H = 0.11
    logg = 4.52

    _, flux, _ = Korg.synth(; Teff=Teff, logg=logg, M_H=M_H, linelist=linelist,
                            wavelengths=synth_wls)
    LSF = Korg.compute_LSF_matrix(synth_wls, obs_wls, R)
    fake_data = LSF * flux

    err = 0.01 * ones(length(obs_wls)) # don't actually apply error to keep tests deterministic

    p0 = (; Teff=4500.0, M_H=-0.2, logg=4.1)

    SUITE["fit"] = BenchmarkGroup()
    SUITE["fit"]["fit_spectrum"] = @benchmarkable Korg.Fit.fit_spectrum($obs_wls, $fake_data,
                                                                        $err, $linelist, $p0;
                                                                        R=$R) setup=(GC.gc()) samples=1
end
