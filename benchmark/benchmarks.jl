using BenchmarkTools, Korg

const SUITE = BenchmarkGroup()

# the PACKAGE_VERSION constant can be used to run tweaked benchmarks for different versions of Korg

# Stellar parameters for benchmarks
const BENCHMARK_PARAMS = [
    ("sun", (Teff=5777, logg=4.44, m_H=0.0)),
    ("M dwarf", (Teff=3500, logg=4.6, m_H=0.0))
]

const WAVELENGTH_RANGES_AND_LINELISTS = [
    ("one wavelength", (5000, 5000), []),
    ("4000 Å - 8000 Å VALD solar", (4000, 8000), Korg.get_VALD_solar_linelist()),
    ("APOGEE DR17 w/o water", (15000, 17000), Korg.get_APOGEE_DR17_linelist(; include_water=false)),
    ("100 Å of APOGEE w/ water", (15500, 15600), Korg.get_APOGEE_DR17_linelist())
]

SUITE["synthesis"] = BenchmarkGroup()
for (i, (param_name, params)) in enumerate(BENCHMARK_PARAMS)
    SUITE["synthesis"][param_name] = BenchmarkGroup()

    # Test different wavelength ranges
    for (j, (range_name, wl_range, linelist)) in enumerate(WAVELENGTH_RANGES_AND_LINELISTS)
        try
            SUITE["synthesis"][param_name][range_name] = @benchmarkable synth(Teff=$(params.Teff),
                                                                              logg=$(params.logg),
                                                                              m_H=$(params.m_H),
                                                                              linelist=$linelist,
                                                                              wavelengths=$wl_range) setup=(GC.gc()) seconds=30
        catch e
            @warn "Could not run benchmark for $param_name/$range_name: $e"
        end
    end
end

results = run(SUITE)
println(results)
