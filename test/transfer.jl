@testset "transfer" begin
    using SpecialFunctions: expint

    @testset "expint" begin
        xs = 2:0.01:8
        @test Korg.RadiativeTransfer.exponential_integral_2.(xs)≈expint.(2, xs) rtol=1e-3
    end

    @testset "compare transfer schemes" begin
        # 3.0 is where we transition from planar to spherical atmospheres in the marcs grid
        atm = interpolate_marcs(5000.0, 3.0)
        patm = Korg.PlanarAtmosphere(atm)

        sol = synthesize(atm, [], format_A_X(), 5000, 5001)
        ref_ind = 1

        wl_cm = sol.wavelengths * 1e-8
        S = Korg.blackbody.(Korg.get_temps(atm), wl_cm')

        for atm in [atm, patm], τ_scheme in ["anchored", "bezier"],
            I_scheme in ["linear", "linear_flux_only", "bezier"],
            include_inward_rays in [true, false]

            flux, _ = Korg.RadiativeTransfer.radiative_transfer(atm, sol.alpha, S, 20;
                                                                include_inward_rays=include_inward_rays,
                                                                α_ref=sol.alpha[:, ref_ind],
                                                                τ_scheme=τ_scheme,
                                                                I_scheme=I_scheme)

            atmtype = if atm isa Korg.PlanarAtmosphere
                "planar"
            else
                "spherical"
            end
            @test assert_allclose_grid(sol.flux, flux,
                                       [("λ [$I_scheme $τ_scheme $atmtype]", sol.wavelengths, "Å")];
                                       rtol=0.05, print_rachet_info=false)
        end
    end

    @testset "mu grid" begin
        @testset "valid specifications" for mu_specification in [
            5,
            20,
            [0, 0.2, 1.0],
            0:0.5:1.0,
            [1.0]
        ]
            # what length should it be?
            if mu_specification isa Number
                n = mu_specification
            else
                n = length(mu_specification)
            end

            μs, ws = Korg.RadiativeTransfer.generate_mu_grid(mu_specification)
            @test length(μs) == n
            @test length(ws) == n
            @test sum(ws) ≈ 1.0
        end

        @testset "invalid mu specifications" begin
            @test_throws ArgumentError Korg.RadiativeTransfer.generate_mu_grid(1:-0.5:0)
            @test_throws ArgumentError Korg.RadiativeTransfer.generate_mu_grid(1:0.5:2)
            @test_throws ArgumentError Korg.RadiativeTransfer.generate_mu_grid([-0.5, 0.5])
        end
    end
end
