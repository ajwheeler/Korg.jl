"""
Compute nₑ (number density of free electrons) in a pure Hydrogen atmosphere, where `nH_tot` is the
total number density of H I and H II (in cm⁻³), the temperature is `T`, and `HI_partition_val` is
the the value of the H I partition function.

This is a relatively naive implementation. More numerically stable solutions exist.
"""
function electron_ndens_Hplasma(nH_tot, T, H_I_partition_val = 2.0)
    # Define the Saha equation as: nₑ*n_{H II} / n_{H I} = RHS
    # coef ∼ 4.829e15
    coef = 2.0 * (2.0*π*Korg.electron_mass_cgs*Korg.kboltz_cgs / Korg.hplanck_cgs^2)^1.5
    RHS = coef * T^1.5 * exp(-Korg.RydbergH_eV/(Korg.kboltz_eV*T))/H_I_partition_val
    # In a pure Hydrogen atmosphere: nₑ = n_{H II}. The Saha eqn becomes:  nₑ²/(nH_tot - ne) = RHS
    # We recast the Saha eqn as: a*nₑ² + b*nₑ + c = 0 and compute the coefficients
    a, b, c = (1.0, RHS, -1*RHS*nH_tot)
    # solve quadratic equation. Since b is always positive and c is always negative:
    #    (-b + sqrt(b²-4*a*c))/(2*a) is always ≥ 0
    #    (-b - sqrt(b²-4*a*c))/(2*a) is always negative
    nₑ = (-b + sqrt(b*b-4*a*c))/(2*a)
    nₑ
end

@testset "stat mech" begin
    @testset "pure Hydrogen atmosphere" begin
        nH_tot = 1e15
        # specify χs and Us to decouple this testset from other parts of the code
        χs = [(Korg.RydbergH_eV, -1.0, -1.0)]
        Us = Dict([Korg.species"H_I"=>(T -> 2.0), Korg.species"H_II"=>(T -> 1.0)])
        # iterate from less than 1% ionized to more than 99% ionized
        for T in [3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 1e4, 1.1e4, 1.2e4, 1.3e4, 1.4e4, 1.5e5]
            nₑ = electron_ndens_Hplasma(nH_tot, T, 2.0)
            wII, wIII = Korg.saha_ion_weights(T, nₑ, 1, χs, Us)
            @test wIII ≈ 0.0 rtol = 1e-15
            rtol = (T == 1.5e5) ? 1e-9 : 1e-14
            @test wII/(1 + wII + wIII) ≈ (nₑ/nH_tot) rtol= rtol
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

    @testset "molecular equilibrium" begin
        #solar abundances
        nX_ntot = @. 10^(Korg.asplund_2009_solar_abundances - 12)
        nX_ntot ./= sum(nX_ntot)

        nₜ = 1e15 
        nₑ_initial_guess = 1e12
        T = 5700
        nₑ, n_dict = Korg.chemical_equilibrium(T, nₜ, nₑ_initial_guess, nX_ntot, 
                                               Korg.ionization_energies, Korg.default_partition_funcs, 
                                               Korg.default_log_equilibrium_constants)

        @test_logs (:warn, r"Electron number density differs") Korg.chemical_equilibrium(
                                            T, nₜ, 1.0, nX_ntot, Korg.ionization_energies, 
                                            Korg.default_partition_funcs, 
                                            Korg.default_log_equilibrium_constants)
        

        # plasma is net-neutral
        positive_charge_density =  mapreduce(+, pairs(n_dict)) do (species, n)
            n * species.charge
        end
        @test nₑ ≈ positive_charge_density rtol=1e-5

        #make sure number densities are sensible
        @test (n_dict[Korg.species"C_III"] < n_dict[Korg.species"C_II"] < n_dict[Korg.species"C_I"] < 
               n_dict[Korg.species"H_II"] < n_dict[Korg.species"H_I"])

        @testset "conservation of nuclei: $(Korg.atomic_symbols[Z])" for Z in 1:Korg.MAX_ATOMIC_NUMBER
            total_n = mapreduce(+, collect(keys(n_dict))) do species
                sum(Korg.get_atoms(species.formula) .== Z) * n_dict[species]
            end
            @test total_n ≈ nX_ntot[Z] * nₜ
        end
    end

    @testset "compare to Barklem and Collet partiion functions" begin
        logTs = 1:0.01:4 #B&C partition functions are only defined up to 10,000 K
        BC_Us = Korg.read_Barklem_Collet_table("data/BarklemCollet2016-atomic_partition.dat");
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

    @testset "ForwardDiff.Dual-specific methods" begin
        A_X = format_A_X()
        abs_abundances = 10 .^ (A_X .- 12)
        abs_abundances ./= sum(abs_abundances)

        g(x) = Korg.chemical_equilibrium(x[1], x[2], x[2]*1e-20, abs_abundances, Korg.ionization_energies, 
                    Korg.default_partition_funcs, Korg.default_log_equilibrium_constants;
                    electron_number_density_warn_threshold=1e100)[1]

        x = [5000 ; 1e12]
        
        g1 = ForwardDiff.gradient(g, x)
        g2 = FiniteDiff.finite_difference_gradient(g, x)
        @test isapprox(g1, g2, rtol=1e-5)

        h1 = ForwardDiff.hessian(g, x)
        h2 = FiniteDiff.finite_difference_hessian(g, x)
        @test isapprox(h1, h2, rtol=1e-5)
    end
end

