@testset "stat mech" begin
    @testset "Saha eqn: pure H" begin
        # Compute nₑ (number density of free electrons) in a pure Hydrogen atmosphere, where `nH_tot` is the
        # total number density of H I and H II (in cm⁻³), the temperature is `T`, and `HI_partition_val` is
        # the the value of the H I partition function.
        # This is a relatively naive implementation. More numerically stable solutions exist.
        function electron_ndens_Hplasma(nH_tot, T, H_I_partition_val=2.0)
            # Define the Saha equation as: nₑ*n_{H II} / n_{H I} = RHS
            # coef ∼ 4.829e15
            coef = 2.0 *
                   (2.0 * π * Korg.electron_mass_cgs * Korg.kboltz_cgs / Korg.hplanck_cgs^2)^1.5
            RHS = coef * T^1.5 * exp(-Korg.RydbergH_eV / (Korg.kboltz_eV * T)) / H_I_partition_val
            # In a pure Hydrogen atmosphere: nₑ = n_{H II}. The Saha eqn becomes:  nₑ²/(nH_tot - ne) = RHS
            # We recast the Saha eqn as: a*nₑ² + b*nₑ + c = 0 and compute the coefficients
            a, b, c = (1.0, RHS, -1 * RHS * nH_tot)
            # solve quadratic equation. Since b is always positive and c is always negative:
            #    (-b + sqrt(b²-4*a*c))/(2*a) is always ≥ 0
            #    (-b - sqrt(b²-4*a*c))/(2*a) is always negative
            nₑ = (-b + sqrt(b * b - 4 * a * c)) / (2 * a)
            nₑ
        end

        nH_tot = 1e15
        # specify χs and Us to decouple this testset from other parts of the code
        χs = [(Korg.RydbergH_eV, -1.0, -1.0)]
        Us = Dict([Korg.species"H_I" => (T -> 2.0), Korg.species"H_II" => (T -> 1.0)])
        # iterate from less than 1% ionized to more than 99% ionized
        for T in [3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 1e4, 1.1e4, 1.2e4, 1.3e4, 1.4e4, 1.5e5]
            nₑ = electron_ndens_Hplasma(nH_tot, T, 2.0)
            wII, wIII = Korg.saha_ion_weights(T, nₑ, 1, χs, Us)
            @test wIII≈0.0 rtol=1e-15
            rtol = (T == 1.5e5) ? 1e-9 : 1e-14
            @test wII/(1+wII+wIII)≈(nₑ/nH_tot) rtol=rtol
        end
    end

    @testset "monotonic N ions Temperature dependence" begin
        weights = [Korg.saha_ion_weights(T, 1.0, 7, Korg.ionization_energies,
                                         Korg.default_partition_funcs) for T in 1:100:10000]
        #N II + NIII grows with T === N I shrinks with T
        @test issorted(first.(weights) + last.(weights))

        # NIII grows with T
        @test issorted(last.(weights))
    end

    @testset "chemical/ionization equilibrium" begin
        #solar abundances
        abs_abundances = @. 10^(Korg.default_solar_abundances - 12)
        abs_abundances ./= sum(abs_abundances)

        # convieniece wrapper to supply default values
        function run_chemical_equilibrium(T, nₜ, nₑ_guess, absolote_abundances)
            Korg.chemical_equilibrium(T, nₜ, nₑ_guess, absolote_abundances,
                                      Korg.ionization_energies,
                                      Korg.default_partition_funcs,
                                      Korg.default_log_equilibrium_constants;
                                      electron_number_density_warn_threshold=Inf)
        end

        @testset "solver is insensitive to the nₑ guess" begin
            # Across many decades of initial guess, all final solutions should match.
            T, nₜ = 3500.0, 1e16
            results = map([1e6, 1e9, 1e11, 1e13, 1e15]) do nₑ_guess
                run_chemical_equilibrium(T, nₜ, nₑ_guess, abs_abundances)
            end
            ref_nₑ, ref_n_dict = results[1]
            for (nₑ, n_dict) in results[2:end]
                @test nₑ≈ref_nₑ rtol=1e-8
                # spot-check a few species across the abundance range
                @test assert_allclose_dict(ref_n_dict, n_dict; rtol=1e-7, print_rachet_info=false)
            end
        end

        @testset "solver convergence is deterministic" begin
            # Running twice from the same inputs should give bitwise identical results
            T, nₜ, nₑ_g = 3000.0, 1e16, 1e10
            nₑ_1, n_dict_1 = run_chemical_equilibrium(T, nₜ, nₑ_g, abs_abundances)
            nₑ_2, n_dict_2 = run_chemical_equilibrium(T, nₜ, nₑ_g, abs_abundances)
            @test nₑ_1 == nₑ_2
            for k in keys(n_dict_1)
                @test n_dict_1[k] == n_dict_2[k]
            end
        end

        @testset "initial nₑ warnings" begin
            T = 5700
            nₜ = 1e15
            nₑ = 1e12
            nₑ_initial_guess = nₜ / 1e3

            args = (T, nₜ, 1.2nₑ, abs_abundances, Korg.ionization_energies,
                    Korg.default_partition_funcs, Korg.default_log_equilibrium_constants)
            @test_logs (:warn, r"Electron number density differs") Korg.chemical_equilibrium(args...)
            # no warning logged if we tweak either the threshold or the minimum applicable nₑ
            @test_logs min_level=Logging.Warn Korg.chemical_equilibrium(args...;
                                                                        electron_number_density_warn_threshold=0.3)
            @test_logs min_level=Logging.Warn Korg.chemical_equilibrium(args...;
                                                                        electron_number_density_warn_min_value=1.01 *
                                                                                                               nₑ)
        end

        # Verifies that chemical_equilibrium results are internally consistent and make sense 
        function check_solution(T, nₑ, n_dict, nₜ, abs_abundances; rtol=1e-8)
            @test isfinite(nₑ) && nₑ > 0
            @test all(v -> isfinite(v) && v >= 0, values(n_dict))

            # charge balance
            positive_charge_density = mapreduce(+, pairs(n_dict)) do (species, n)
                n * species.charge
            end
            @test nₑ≈positive_charge_density rtol=rtol

            # spot-check H-
            n_Hminus = n_dict[Korg.species"H-"]
            n_HI = n_dict[Korg.species"H_I"]
            @test n_Hminus≈Korg.Hminus_nK(T)*n_HI*nₑ rtol=1e-5
            @test n_Hminus > 0

            # total number density is correct
            total_number_density = nₑ + sum(values(n_dict))
            @test total_number_density≈nₜ rtol=rtol

            # abundances in == abundances out
            nuclear_n_densities = map(1:Korg.MAX_ATOMIC_NUMBER) do Z
                mapreduce(+, collect(keys(n_dict))) do species
                    sum(Korg.get_atoms(species.formula) .== Z) * n_dict[species]
                end
            end
            nuclear_n_densities ./= sum(nuclear_n_densities)
            @test nuclear_n_densities≈abs_abundances rtol=rtol

            # n(C_I) < n(H_I) always holds (carbon ≪ hydrogen by abundance regardless of ionization).
            @test n_dict[Korg.species"C_I"] < n_dict[Korg.species"H_I"]

            # if molecules aren't in play, n(C I) < n(C II) < n(C III) because of Saha
            if T ≤ 6000
                @test n_dict[Korg.species"C_III"] < n_dict[Korg.species"C_II"] <
                      n_dict[Korg.species"C_I"]
            end
        end

        @testset "all layers of representative atmospheres" begin
            cases = [
                # (Teff, logg, [M/H], [α/M], [C/M],  label)
                (5777.0, 4.44, 0.0, 0.0, 0.0, "solar"),
                (4500.0, 2.5, -0.5, 0.2, 0.0, "K giant"),
                (3500.0, 4.8, 0.0, 0.0, 0.0, "M dwarf"),
                (2800.0, 0.5, -2.5, -1.0, -1.5, "extreme cool C-poor giant"),
                (2800.0, 4.5, -2.5, 0.0, 0.0, "cool metal-poor dwarf"),
                (4000.0, 0.5, 0.5, 0.0, 0.0, "cool metal-rich giant"),
                (3000.0, 5.0, -4.0, 0.4, 0.0, "low-Z cool dwarf"),
                (2500.0, 5.0, -5.0, 0.4, 0.0, "low-Z extreme cool dense dwarf")
            ]
            @testset "$(lbl)" for (Teff, logg, M_H, α_M, C_M, lbl) in cases
                A_X = Korg.format_A_X(M_H, α_M + M_H, Dict("C" => C_M + M_H))
                atm = Korg.interpolate_marcs(Teff, logg, A_X; clamp_abundances=true)

                ab = 10 .^ (A_X .- 12)
                ab ./= sum(ab)
                for layer in atm.layers
                    T = layer.temp
                    nₜ = layer.number_density
                    nₑ_atm = layer.electron_number_density
                    nₑ, n_dict = run_chemical_equilibrium(T, nₜ, nₑ_atm, ab)
                    check_solution(T, nₑ, n_dict, nₜ, ab)
                end
            end
        end

        @testset "C/O ratio sweep" begin
            # When C/O crosses 1, the dominant carbon and oxygen carriers flip quickly, 
            # which is hard for the chemical equilibrium solver.
            T, nₜ, nₑg = 1500.0, 1e17, 1e9
            A_X0 = Korg.format_A_X()
            # vary C abundance, keep O fixed
            log_CO_ratios = [-0.5, -0.2, -0.1, -0.02, 0.0, 0.02, 0.1, 0.2, 0.5]
            for Δ in log_CO_ratios
                A_X = copy(A_X0)
                A_X[6] = A_X[8] + Δ   # log10(C/O) = Δ
                ab = 10 .^ (A_X .- 12)
                ab ./= sum(ab)
                nₑ, n_dict = run_chemical_equilibrium(T, nₜ, nₑg, ab)
                check_solution(T, nₑ, n_dict, nₜ, ab)

                # CO can't exceed the limiting element's total nuclear density (sum atoms-in-species
                # over all species containing C or O, respectively).
                n_C_total = sum(n_dict[s] * count(==(6), Korg.get_atoms(s.formula))
                                for s in keys(n_dict))
                n_O_total = sum(n_dict[s] * count(==(8), Korg.get_atoms(s.formula))
                                for s in keys(n_dict))
                @test n_dict[Korg.species"CO"] <= min(n_C_total, n_O_total)
            end
        end

        @testset "ForwardDiff.Dual-specific methods" begin
            g(x) = run_chemical_equilibrium(x[1], x[2], x[2] * 1e-20, abs_abundances)[1]

            for x in ([5000.0, 1e15], [3000.0, 1e16], [2500.0, 1e17])
                g1 = ForwardDiff.gradient(g, x)
                g2 = FiniteDiff.finite_difference_gradient(g, x)
                @test isapprox(g1, g2, rtol=1e-5)

                h1 = ForwardDiff.hessian(g, x)
                h2 = FiniteDiff.finite_difference_hessian(g, x)
                @test isapprox(h1, h2, rtol=1e-5)
            end
        end
    end

    @testset "compare to Barklem and Collet partiion functions" begin
        logTs = 1:0.01:4 #B&C partition functions are only defined up to 10,000 K
        BC_Us = Korg.read_Barklem_Collet_table("data/BarklemCollet2016-atomic_partition.dat")
        Ts = 10 .^ logTs

        # I'm only comparing the first 10 neutal species because of unexplained weirdness with the
        # B&C numbers.
        @testset for spec in Korg.Species.(Korg.atomic_symbols[1:10] .* " I")
            BC_U = BC_Us[spec].(logTs * log(10))
            korg_U = Korg.default_partition_funcs[spec].(logTs * log(10))

            # we can't get closer than 2% because NIST has changed some things
            # (There may also be edge-case energy levels that are weird in some way that are only
            # included in one calculation, but I think this effect is small.)
            @test assert_allclose_grid(korg_U, BC_U, [("T", Ts, "K")]; rtol=0.01,
                                       print_rachet_info=false)
        end
    end

    @testset "atomic data" begin
        @test (Korg.MAX_ATOMIC_NUMBER
               == length(Korg.atomic_masses)
               == length(Korg.asplund_2009_solar_abundances)
               == length(Korg.asplund_2020_solar_abundances)
               == length(Korg.grevesse_2007_solar_abundances)
               == length(Korg.bergemann_2025_solar_abundances)
               == length(Korg.magg_2022_solar_abundances))

        @test (Korg.get_mass(Korg.Formula("CO")) ≈
               Korg.get_mass(Korg.Formula("C")) + Korg.get_mass(Korg.Formula("O")))
        @test Korg.get_mass(Korg.Formula("C2")) ≈ 2Korg.get_mass(Korg.Formula("C"))
    end

    @testset "ionization energies" begin
        @test length(Korg.ionization_energies) == 92
        @test Korg.ionization_energies[Korg.atomic_numbers["H"]] == [13.5984, -1.000, -1.000]
        @test Korg.ionization_energies[Korg.atomic_numbers["Ru"]] == [7.3605, 16.760, 28.470]
        @test Korg.ionization_energies[Korg.atomic_numbers["U"]] == [6.1940, 11.590, 19.800]
    end

    @testset "O I-III and CN partition functions are (nearly) monotonic in T" begin
        lnTs = 0:0.1:log(100_000.0)
        nearly_monotonic(Us) = all(diff(Us) .> -1e-4)
        @test nearly_monotonic(Korg.default_partition_funcs[Korg.species"O_I"].(lnTs))
        @test nearly_monotonic(Korg.default_partition_funcs[Korg.species"O_II"].(lnTs))
        @test nearly_monotonic(Korg.default_partition_funcs[Korg.species"O_III"].(lnTs))
        @test nearly_monotonic(Korg.default_partition_funcs[Korg.species"CN_I"].(lnTs))
    end
end
