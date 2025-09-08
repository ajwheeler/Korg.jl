@testset "synth" begin
    default_linelist = Korg.get_VALD_solar_linelist()
    default_atm = interpolate_marcs(5000, 4.5, format_A_X())
    default_ps = (; Teff=5000, logg=4.5)

    @testset "basic functionality" begin
        # Compare with direct Korg.synthesize call
        wls, flux, cntm = synth(; default_ps..., wavelengths=(5000, 5001))
        sol = synthesize(default_atm, default_linelist, format_A_X(), (5000, 5001))
        @test wls == sol.wavelengths
        @test flux ≈ sol.flux ./ sol.cntm
        @test cntm ≈ sol.cntm
    end

    @testset "rectification" begin
        wls1, flux_rect, cntm = synth(; default_ps..., wavelengths=(5000, 5001), rectify=true)
        wls2, flux_raw, _ = synth(; default_ps..., wavelengths=(5000, 5001), rectify=false)
        @test flux_rect ≈ flux_raw ./ cntm
    end

    @testset "resolution" begin
        wls, flux_lsf, _ = synth(; default_ps..., wavelengths=(5000, 5001), R=50_000)
        wls_raw, flux_raw, _ = synth(; default_ps..., wavelengths=(5000, 5001), R=Inf)
        @assert wls == wls_raw
        flux_manual_lsf = Korg.apply_LSF(flux_raw, wls_raw, 50_000)
        @test flux_lsf ≈ flux_manual_lsf
    end

    @testset "rotation" begin
        wl_spec = [(5000, 5001), (5002, 5003)]
        wls, flux_rot, _ = synth(; default_ps..., wavelengths=wl_spec, vsini=10)
        wls_raw, flux_raw, _ = synth(; default_ps..., wavelengths=wl_spec, vsini=0)
        @assert wls == wls_raw
        flux_manual_rot = Korg.apply_rotation(flux_raw, wl_spec, 10)
        @test flux_rot ≈ flux_manual_rot
    end

    @testset "abundances" begin
        # Test that abundance changes match direct synthesize calls
        wls, flux, _ = synth(; default_ps..., wavelengths=(5000, 5001), M_H=-0.5, alpha_H=0.4,
                             C=-0.3)

        atm = interpolate_marcs(5000, 4.5, format_A_X(-0.5, 0.4, Dict("C" => -0.3)))
        sol = synthesize(atm, default_linelist, format_A_X(-0.5, 0.4, Dict("C" => -0.3)),
                         (5000,
                          5001))

        @test wls == sol.wavelengths
        @test flux ≈ sol.flux ./ sol.cntm
    end

    @testset "atmospheric parameters" begin
        # Test that Teff and logg changes match direct synthesize calls
        wls, flux, _ = synth(; default_ps..., wavelengths=(5000, 5001), Teff=4500, logg=2.5)

        atm = interpolate_marcs(4500, 2.5, format_A_X())
        sol = synthesize(atm, default_linelist, format_A_X(), (5000, 5001))

        @test wls == sol.wavelengths
        @test flux ≈ sol.flux ./ sol.cntm
    end

    @testset "microturbulence" begin
        wls, flux, _ = synth(; default_ps..., wavelengths=(5000, 5001), vmic=2.0)
        sol = synthesize(default_atm, default_linelist, format_A_X(), (5000, 5001); vmic=2.0)
        @test flux ≈ sol.flux ./ sol.cntm
    end

    @testset "linelist" begin
        linelist = []
        wls, flux, _ = synth(; default_ps..., wavelengths=(5000, 5001), linelist=linelist)
        sol = synthesize(default_atm, linelist, format_A_X(), (5000, 5001))
        @test flux ≈ sol.flux ./ sol.cntm
    end

    @testset "parameter passing" begin
        # Test that synthesize_kwargs are properly passed through
        wls, flux, _ = synth(; default_ps..., wavelengths=(5000, 5001),
                             synthesize_kwargs=Dict(:vmic => 2.0))

        sol = synthesize(default_atm, default_linelist, format_A_X(), (5000, 5001); vmic=2.0)

        @test wls == sol.wavelengths
        @test flux ≈ sol.flux ./ sol.cntm
    end

    @testset "m_H deprecation error" begin
        msg = "m_H is no longer a supported keyword argument"
        @test_throws msg synth(; default_ps..., wavelengths=(5000, 5001), m_H=0.0)
    end

    @testset "invalid arguments" begin
        msg = "invalid_arg was passed as a keyword argument to synth, but it is not a valid keyword argument."
        @test_throws msg synth(; default_ps..., wavelengths=(5000, 5001), invalid_arg=0.0)

        # be helpful if user is trying to specify an elemental abundance but using the wrong notation
        msg = "Ca instead of Ca_H"
        @test_throws msg synth(; default_ps..., wavelengths=(5000, 5001), Ca_H=0.0)
    end
end
