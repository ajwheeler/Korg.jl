
# to test the implementation of free-free linear absorption coefficients from John75a and John75b,
# we compare the data they provide for He⁻ ff absorption, against our other He⁻ free-free data
function _john75_Heminus_ff(ν_arr, T, nHe_I, ne)
    map(ν_arr) do ν
        Korg.ContinuumAbsorption._combined_john_neg_ion_absorption(
            ν, T, nHe_I, ne, Korg.ContinuumAbsorption._short_wavelength_Heminus_ff_interp,
            Korg.ContinuumAbsorption._long_wavelength_ff_interps["He"]
        )
    end
end


@testset "neutral ff absorption" begin
    @testset "absorption with data from John 75a and John 75b" begin
        # compare absorption coefficients between Korg.ContinuumAbsorption.Heminus_ff and what we
        # would get if we used the data from John 1975 a and John 1975 b. 

        # this spans over the interval for which Korg.ContinuumAbsorption.Heminus_ff is valid
        λ_arr_μm = [0.5063, 1.0, 5.0, 10.0, 11.0, 15.1878]
        ν_arr = Korg.c_cgs ./ (λ_arr_μm ./ 1e4)

        T_val = 3000.0 # Kelvin

        nₑ = 1e13 # cm⁻³
        nHe_I = 1e13 # cm⁻³

        ref = Korg.ContinuumAbsorption.Heminus_ff(ν_arr, T_val, nHe_I, nₑ)
        @test assert_allclose(_john75_Heminus_ff(ν_arr, T_val, nHe_I, nₑ), ref;
                              rtol = 0.2, atol =0.0);
    end

    @testset "N⁻ ff absorption test against data from John 75a" begin
        # to check that Nminus_ff is implemented correctly, let's compare the absorption
        # coefficient against the coefficient computed from John 75a for long wavelengths

        # problem parameters
        λ_arr_μm = [10.0, 11.0, 15.1878]
        T_val = 3000.0
        nₑ = 1e13 # cm⁻³
        nN_I = 1e13 # cm⁻³

        # compute the "reference α" based on John 75a (includes correction for stimulated emission)
        λ_arr_cm = λ_arr_μm ./ 1e4
        Pₑ = nₑ * Korg.kboltz_cgs * T_val
        interp = Korg.ContinuumAbsorption._long_wavelength_ff_interps["N"]
        ref_α = (λ_arr_cm .* λ_arr_cm .* 1e16) .* (interp(T_val) * Pₑ * nN_I)

        # compute the "actual α"
        calc_α = Korg.ContinuumAbsorption.Nminus_ff(Korg.c_cgs ./ (λ_arr_μm ./ 1e4), T_val, nN_I,
                                                    nₑ)

        # there's a factor of ~8 difference in the results. This seems a little extreme.
        # Admittedly, John 1975 speculates that the error on their values at 10 μm might be a
        # factor of 2. But, this still seems disturbingly big...
        #
        # A better test might use reference data from John & Williams 1977
        @test assert_allclose(calc_α, ref_α; rtol = 6.8, atol =0.0);
    end
end
