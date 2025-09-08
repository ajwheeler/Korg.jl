let # for namespace hygiene
    # Stellar parameters for benchmarks
    benchmark_params = [
        ("sun", (Teff=5777, logg=4.44, M_H=0.0)),
        ("M dwarf", (Teff=3500, logg=4.6, M_H=0.0))
    ]

    wavelength_ranges_and_linelists = [
        ("one wavelength", (5000, 5000), []),
        ("4000 Å - 8000 Å VALD solar", (4000, 8000), Korg.get_VALD_solar_linelist()),
        ("APOGEE DR17 w/o water", (15000, 17000),
         Korg.get_APOGEE_DR17_linelist(; include_water=false)),
        ("100 Å of APOGEE w/ water", (15500, 15600), Korg.get_APOGEE_DR17_linelist())
    ]

    SUITE["synthesis"] = BenchmarkGroup()
    for (i, (param_name, params)) in enumerate(benchmark_params)
        SUITE["synthesis"][param_name] = BenchmarkGroup()

        # Test different wavelength ranges
        for (j, (range_name, wl_range, linelist)) in enumerate(wavelength_ranges_and_linelists)
            try
                SUITE["synthesis"][param_name][range_name] = @benchmarkable synth(;
                                                                                  Teff=$(params.Teff),
                                                                                  logg=$(params.logg),
                                                                                  M_H=$(params.M_H),
                                                                                  linelist=$linelist,
                                                                                  wavelengths=$wl_range) setup=(GC.gc()) seconds=30 samples=1
            catch e
                @warn "Could not run benchmark for $param_name/$range_name: $e"
            end
        end
    end
end
