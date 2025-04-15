@testset "autodiff" begin
    using ForwardDiff, FiniteDiff

    linelist = [Korg.get_VALD_solar_linelist(); Korg.get_APOGEE_DR17_linelist()]
    #cover the optical and the IR to catch different H lines
    wls = [6564:0.01:6565, 15_045:0.01:15_046]
    for atm in read_model_atmosphere.([
                                          "data/sun.mod",
                                          "data/s6000_g+1.0_m0.5_t05_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod"
                                      ])
        # the second model atmosphere happens to be in a weird (probably unphysical) part of
        # parameter space where the electron number densities calculated doesn't match the marcs
        # numbers.
        function flux(p)
            synthesize(atm, linelist, format_A_X(p[1], Dict("Ni" => p[2])), wls;
                       vmic=p[3], electron_number_density_warn_threshold=Inf).flux
        end
        #make sure this works.
        J = ForwardDiff.jacobian(flux, [0.0, 0.0, 1.5])
        @test .!any(isnan.(J))
    end

    @testset "autodiff just one abundance" begin
        atm = read_model_atmosphere("data/sun.mod")
        linelist = [Korg.Line(6000e-8, 0.0, Korg.species"C I", 0.0)]
        # If this line is super weak (or removed), the test will fail due to numerics in the
        # abundances and continuum opacities. It used to be a fake line at 5000, but the fact that
        # that's the reference wavelength caused instability in the finie differnces calculation.
        function f(A_C)
            A_X = format_A_X(Dict("C" => A_C))
            synthesize(atm, linelist, A_X, 6000, 6000).flux[1]
        end
        @test FiniteDiff.finite_difference_derivative(f, 0.0)≈ForwardDiff.derivative(f, 0.0) rtol=1e-4
    end

    @testset "line params" begin
        atm = read_model_atmosphere("data/sun.mod")
        function f(loggf)
            linelist = [Korg.Line(6000e-8, loggf, Korg.species"Na I", 0.0)]
            # This used to be a fake line at 5000, but the fact that that's the reference wavelength
            # caused instability in the finie differnces calculation.
            synthesize(atm, linelist, format_A_X(), 6000, 6000).flux[1]
        end
        @test FiniteDiff.finite_difference_derivative(f, 0.0)≈ForwardDiff.derivative(f, 0.0) rtol=1e-4
    end
end
