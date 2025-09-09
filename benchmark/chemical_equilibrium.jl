let
    SUITE["chemical_equilibrium"] = BenchmarkGroup()

    nX_ntot = @. 10^(Korg.default_solar_abundances - 12)
    nX_ntot ./= sum(nX_ntot)

    # T, nₜ, nₑ
    params = [
        (7875.9, 1.40534e17, 6.68689e14),
        (4549.37, 6.92348e15, 5.47132e11),
        (2222.51, 3.95535e14, 5.15098e8)
    ]

    for (T, nₜ, nₑ) in params
        SUITE["chemical_equilibrium"]["chemical_equilibrium ($T K)"] = @benchmarkable Korg.chemical_equilibrium($T,
                                                                                                                $nₜ,
                                                                                                                $nₑ,
                                                                                                                $nX_ntot,
                                                                                                                $Korg.ionization_energies,
                                                                                                                $Korg.default_partition_funcs,
                                                                                                                $Korg.default_log_equilibrium_constants) setup=(GC.gc()) samples=10
    end
end
