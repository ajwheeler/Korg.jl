@testset "atmosphere" begin
    function assert_atmospheres_close(atm1, atm2; rtol=1e-3)
        val = true
        for f in [Korg.get_tau_5000s, Korg.get_zs, Korg.get_temps, Korg.get_number_densities, Korg.get_electron_number_densities]
            val &= assert_allclose(f(atm1), f(atm2); rtol=rtol, print_rachet_info=false)
        end
        val
    end

    @testset "plane-parallel atmosphere" begin
        #the MARCS solar model atmosphere
        atm = Korg.read_model_atmosphere("data/sun.mod")
        @test atm isa Korg.PlanarAtmosphere
        @test length(atm.layers) == 56
        @test issorted([l.temp for l in atm.layers])
        @test atm.layers[1].tau_5000 ≈ 0.00001209483645
        @test atm.layers[1].z == 6.931E+07
        @test atm.layers[1].temp == 4066.8
        @test atm.layers[1].electron_number_density ≈ 3.769664452210607e10
        @test atm.layers[1].number_density ≈ 4.75509171357701e14

        # just make sure these don't error
        Korg.get_tau_5000s(atm)
        Korg.get_zs(atm)
        Korg.get_temps(atm)
        Korg.get_electron_number_densities(atm)
        Korg.get_number_densities(atm)
    end
    @testset "spherical atmosphere" begin
        atm = Korg.read_model_atmosphere(
                "data/s6000_g+1.0_m0.5_t05_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
        @test atm isa Korg.ShellAtmosphere
        @test length(atm.layers) == 56 
        @test issorted([l.temp for l in atm.layers])
        @test atm.R == 2.5827E+12
        @test atm.layers[1].tau_5000 ≈ 4.584584692493259e-5
        @test atm.layers[1].z == 2.222e11
        @test atm.layers[1].temp == 3935.2
        @test atm.layers[1].electron_number_density ≈ 1.7336231777439526e8
        @test atm.layers[1].number_density ≈ 1.5411190391302566e12
    end

    @testset "atmosphere type conversion" begin
        atm = Korg.read_model_atmosphere("data/sun.mod")
        atm2 = Korg.PlanarAtmosphere(Korg.ShellAtmosphere(atm, 7e10)) #arbitrary radius
        @test atm.layers == atm2.layers

        atm = Korg.read_model_atmosphere(
                "data/s6000_g+1.0_m0.5_t05_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
        atm2 = Korg.ShellAtmosphere(Korg.PlanarAtmosphere(atm), 1.0)
        @test [l.tau_5000 for l in atm.layers]                == [l.tau_5000 for l in atm2.layers]
        @test [l.z for l in atm.layers]                       == [l.z for l in atm2.layers]
        @test [l.temp for l in atm.layers]                    == [l.temp for l in atm2.layers]
        @test [l.number_density for l in atm.layers]          == [l.number_density for l in atm2.layers]
        @test [l.electron_number_density for l in atm.layers] == [l.electron_number_density for l in atm2.layers]
    end


    @testset "model atmosphere interpolation" begin
        @testset "methods are equivalent" begin
            teff = 5000.0
            logg = 4.0
            m_H = 0.1
            alpha_H = 0.2
            C_H = 0.3
            atm1 = interpolate_marcs(teff, logg, m_H, alpha_H - m_H, C_H - m_H)
            A_X = format_A_X(m_H, alpha_H, Dict("C"=>C_H); solar_abundances=Korg.grevesse_2007_solar_abundances)
            atm2 = interpolate_marcs(teff, logg, A_X)

            @test assert_atmospheres_close(atm1, atm2; rtol=1e-5)
        end

        @testset "clamping abundances" begin
            m_H_nodes = Korg._sdss_marcs_atmospheres[1][3]

            atm1 = interpolate_marcs(5000.0, 3.0, m_H_nodes[1])
            A_X_2 = format_A_X(m_H_nodes[1] - 1; solar_abundances=Korg.grevesse_2007_solar_abundances)
            atm2 = interpolate_marcs(5000.0, 3.0, A_X_2; clamp_abundances=true)
            @test assert_atmospheres_close(atm1, atm2; rtol=1e-2)

            atm1 = interpolate_marcs(5000.0, 3.0, m_H_nodes[end])
            A_X_2 = format_A_X(m_H_nodes[end] + 1; solar_abundances=Korg.grevesse_2007_solar_abundances)
            atm2 = interpolate_marcs(5000.0, 3.0, A_X_2; clamp_abundances=true)
            @test assert_atmospheres_close(atm1, atm2; rtol=1e-2)
        end

        @testset "grid points" begin
            using Interpolations: linear_interpolation
            # calling the interpolator on grid points should return the same atmosphere
            atm1 = Korg.read_model_atmosphere("data/s5000_g+3.0_m1.0_t02_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
            atm2 = Korg.interpolate_marcs(5000, 3.0)

            # values are not precisely identical.  I think this is due to slightly different solar mixtures.
            @test assert_atmospheres_close(atm1, atm2; rtol=2e-3)

            # two interpolation schemes should give the same result at grid points
            atm1 = Korg.interpolate_marcs(3000, 4.0; resampled_cubic_for_cool_dwarfs=true)
            atm2 = Korg.interpolate_marcs(3000, 4.0; resampled_cubic_for_cool_dwarfs=false)

            logτs1 = log10.(Korg.get_tau_5000s(atm1))
            logτs2 = log10.(Korg.get_tau_5000s(atm2))

            for f in [Korg.get_temps, Korg.get_number_densities, Korg.get_electron_number_densities]
                itp = linear_interpolation(logτs1, f(atm1))
                fs = itp(logτs2[5:end-5])
                @test assert_allclose(f(atm2)[5:end-5], fs; rtol=1e-2)
            end
        end
    end
end
