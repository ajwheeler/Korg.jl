using Interpolations: bounds

@testset "continuum absorption" begin
    #helper stuff
    include("absorption_comparison_funcs.jl")

    #test Peach ff absorption
    include("ff.jl")

    function _bell_berrington_87_Hminus_ff_absorption_coef(ν, T)
        # Compute the quantity that Bell & Berrington (1987) refer to as the "absorption coefficient"
        # in units of 1e-26 cm⁴ dyn⁻¹ (so that we can compare against their table). Essentially, we are
        # just solving for:
        #     ff_linear_abs_coef(H⁻) / (n(H I) * Pₑ)
        # This could be viewed as the free-free cross-section (including the correction for stimulated
        # stimulated emission) per partial electron pressure.

        # The choices of nH_I, ne, and ρ are completely unimportant since α_ff(H⁻) has no explicit
        # dependence on any of them.
        nH_I, ne = 1e17, 1e14

        # To allow the best comparison to real data let's set the H I partition function to 2, so that
        # n(H I, n = 1) = n(H I).
        nH_I_div_partition = nH_I / 2.0

        linear_abs_coef = Korg.ContinuumAbsorption.Hminus_ff([ν], T, nH_I_div_partition, ne)[1]
        Pₑ = ne * Korg.kboltz_cgs * T
        linear_abs_coef / (nH_I * Pₑ * 1e-26)
    end

    # the first axis of ref_val_matrix should have the same length as λ_vals_cm
    # the second axis of ref_val_matrix should have the same length as T_vals
    function _compare_against_table(ref_val_matrix, λ_vals_cm, T_vals, calc_func, target_precision,
                                    verbose=true)
        @assert ((length(size(ref_val_matrix)) == 2) &
                 (length(size(λ_vals_cm)) == 1) &
                 (length(size(T_vals)) == 1))
        @assert size(ref_val_matrix)[1] == size(λ_vals_cm)[1]
        @assert size(ref_val_matrix)[2] == size(T_vals)[1]

        ν_vals = Korg.c_cgs ./ λ_vals_cm

        success = true
        for index in 1:length(T_vals)
            cur_T = T_vals[index]
            ref_vals = view(ref_val_matrix, :, index)
            calc_vals = calc_func.(ν_vals, cur_T)
            precision = (abs.(calc_vals .- ref_vals) ./ ref_vals)
            max_err, max_err_ind = findmax(precision)
            if max_err > target_precision
                if verbose
                    if success
                        @info string("There is a problem. The max error of should be: ",
                                     target_precision * 100, "%.")
                    end
                    @info string("Max error for T = ", cur_T, "K is: ", 100 * max_err,
                                 "% at λ = ", (λ_vals_cm[max_err_ind] * 1.e8), " Å.")
                end
                success = false
            end
        end
        return success
    end

    """
    This is taken from equation 8.13 of Gray (2005), which is a a polynomial fit against Table 1
    of Bell & Berrington (1987) [https://doi.org/10.1088/0022-3700/20/4/019] (the source Korg uses).
    The polynomial fits `f(H⁻)` in the range 2520 K ≤ T ≤ 10080 K and 2604 Å ≤ λ ≤ 113918 Å.
    According to Grey, the polynomial fit to the tabulated data typically has 1% precision. We find
    that at worst, the discrepancy never exceeds 2.25%.
    """
    function _Hminus_ff_cross_section_Gray(λ, logθ, nₑ)
        logλ = log10(λ)
        log2λ = logλ * logλ
        log3λ = log2λ * logλ
        log4λ = log3λ * logλ

        f0 = -2.2763 - 1.6850 * logλ + 0.76661 * log2λ - 0.053346 * log3λ
        f1 = +15.2827 - 9.2846 * logλ + 1.99381 * log2λ - 0.142631 * log3λ
        f2 = -197.789 + 190.266 * logλ - 67.9775 * log2λ + 10.6913 * log3λ - 0.625151 * log4λ

        Pe = Korg.kboltz_cgs * (5040 ./ (10^logθ)) * nₑ
        1e-26 * Pe * 10.0^(f0 + f1 * logθ + f2 * (logθ * logθ))
    end

    @testset "H⁻ free-free absorption" begin
        @testset "compare to Gray polynomial" begin
            # these are from the  Bell & berrrington table, though we could use any reasonable vals
            # we cut off vals outside support of polynomial approximation
            # each row corresponds to a separate λ (in Å)
            # each column corresponds to a separate θ value where θ = 5040/T
            λ_vals = [9113, 7595, 6510, 5696, 5063, 4557, 3645, 3038, 2604]
            ν_vals = Korg.c_cgs ./ (λ_vals * 1e-8)
            θ_vals = [0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8]
            T_vals = 5040.0 ./ θ_vals

            #assume nₑ = n(H I) = 1.0
            α_ref = _Hminus_ff_cross_section_Gray.(λ_vals, log10.(θ_vals)', 1.0)
            #use the unwrapped (not bounds-checked) version for easier broadcasting
            n_H_I_div_U = 0.5
            α_korg = Korg.ContinuumAbsorption._Hminus_ff.(ν_vals, T_vals', n_H_I_div_U, 1.0)

            @test assert_allclose_grid(α_korg, α_ref, [("λ", λ_vals, "Å"), ("T", T_vals, "K")];
                                       rtol=0.0225)
        end

        # we only have measurements for b and c
        @testset "Gray (2005) Fig 8.5$panel comparison" for panel in ["b", "c"]
            calculated, ref = Gray_opac_compare.Gray05_comparison_vals(panel, "Hminus_ff")
            @test all(calculated .≥ 0.0)
            @test assert_allclose(calculated, ref; rtol=0, atol=0.04)
        end
    end

    """
        _Hminus_bf_cross_section_Gray(λ)

    Compute the H⁻ bound-free cross-section at a given wavelength (specified in Å). The cross-section
    is converted to units of cm^2 per H⁻ particle and does NOT include a correction for stimulated
    emission.

    This function uses the polynomial provided in equation 8.11 of Gray (2005), that fits the tabulated
    data from Wishart (1979). While Gray (2005) claims that the polynomial fits the data with 0.2%
    precision for 2250 Å ≤ λ ≤ 15000 Å, in practice we find that it fits the data to better than 0.25%
    precision.  It fits the McLaughlin (2017) data (what Korg uses) within 3%.
    """
    function _Hminus_bf_cross_section_Gray(λ::Real)
        # we need to somehow factor out this bounds checking
        if !(2250 <= λ <= 15000.0)
            throw(DomainError(λ, "The wavelength must lie in the interval [2250 Å, 15000 Å]"))
        end

        λ2 = λ * λ
        λ3 = λ * λ2
        λ4 = λ * λ3
        λ5 = λ * λ4
        λ6 = λ * λ5

        1e-18 * (1.99654 - 1.18267e-5 * λ + 2.64243e-6 * λ2 - 4.40524e-10 * λ3 + 3.23992e-14 * λ4
         -
         1.39568e-18 * λ5 + 2.78701e-23 * λ6)
    end

    @testset "H⁻ bound-free absorption" begin
        @testset "comparison to Gray" begin
            # compare the tabulated data from McLaughlin (2017) to the polynomial fits from Gray (2005)
            λs = 2250:1.0:15000
            νs = Korg.c_cgs ./ (λs * 1e-8)
            grayvals = _Hminus_bf_cross_section_Gray.(λs)
            korgvals = Korg.ContinuumAbsorption._Hminus_bf_cross_section.(νs)
            @test assert_allclose_grid(korgvals, grayvals, [("λ", λs, "Å")], rtol=0.02)
        end

        @testset "Gray (2005) Fig 8.5$panel comparison" for panel in ["a", "b", "c"]
            calculated, ref = Gray_opac_compare.Gray05_comparison_vals(panel, "Hminus_bf")
            @test all(calculated .≥ 0.0)
            @test assert_allclose(calculated, ref; rtol=0.03)
        end

        @testset "extreme values" begin
            #test frequencies outside of the range where tabulated values are available

            #should return 0 when below the electron affinity (ionization energy)
            ν_small = Korg.c_cgs / (20_000 * 1e-8) #20,000 Å
            @test Korg.ContinuumAbsorption._Hminus_bf_cross_section(ν_small) == 0

            #cross section should monotonically increase with increasing ν / decreasing λ as we go from
            #the lowest energy with a non-0 cross section, through the region without tabulated values,
            #to the lowest tabulated value
            νs = ((Korg.ContinuumAbsorption._H⁻_ion_energy/Korg.hplanck_eV)
                  :1e12:
                  (Korg.ContinuumAbsorption._min_H⁻_interp_ν+2e13))
            σs = Korg.ContinuumAbsorption._Hminus_bf_cross_section.(νs)
            @test all(diff(σs) .> 0)
        end
    end

    """
    This follows equation 8.16 of Grey (2005) which provides a polynomial fit absorption to data
    tabulated in John (1994) [https://ui.adsabs.harvard.edu/abs/1994MNRAS.269..871J/abstract] (the
    source used by Korg).

    According to that equation, the "continuous absorption coefficient per Hydrogen" is given by:
    f(He⁻_ff)*Pₑ*A(He) / (1 + Φ(He)/Pₑ)
    in which

      - log_10(f(He⁻_ff)) is the polynomial term (Grey actually uses the variable α in place of f)
      - A(He) is the number abundance of Helium particles relative to Hydrogen particles. For reference,
        A(He) = [n(He I) + n(He II) + n(He III)] / [n(H I) + n(H II)]
      - Pₑ is the partial pressure contributed by electrons, Pₑ = nₑkT.
      - `1/(1 + Φ(He)/Pₑ)` comes from the Saha equation and expresses `n(He I) / [n(He I) + n(He II)]`.
        In the above expression, f(He⁻_ff)*Pₑ specifies the free-free atomic absorption coefficient per
        ground state He I atom.  This is what we return.

    For 5063 Å ≤ λ ≤ 151878 Å and 2520 K ≤ T ≤ 10080 K, Gray (2005) claims that the polynomial fit the
    tabulated data at a precision of 1% or better. In practice, we found that it only fits the data to
    better than 3.1% (it's possible that for smaller λ the fit may be better). For reference, the
    tabulated data in these ranges of values consist of an irregularly spaced rectangular grid with 15
    λ values and 9 Temperature values. According to John (1994), improved calculations are unlikely to
    alter the tabulated data for λ > 1e4Å, "by more than about 2%." The errors introduced by the
    approximations for 5.06e3 Å ≤ λ ≤ 1e4 Å "are expected to be well below 10%."
    """
    function _Heminus_ff_cross_section_Gray(λ, θ, nₑ)
        θ2 = θ * θ
        θ3 = θ2 * θ
        θ4 = θ3 * θ

        c0 = 9.66736 - 71.76242 * θ + 105.29576 * θ2 - 56.49259 * θ3 + 10.69206 * θ4
        c1 = -10.50614 + 48.28802 * θ - 70.43363 * θ2 + 37.80099 * θ3 - 7.15445 * θ4
        c2 = 2.74020 - 10.62144 * θ + 15.50518 * θ2 - 8.33845 * θ3 + 1.57960 * θ4
        c3 = -0.19923 + 0.77485 * θ - 1.13200 * θ2 + 0.60994 * θ3 - 0.11564 * θ4

        logλ = log10(λ)

        # f includes contribution from stimulated emission
        f = 1e-26 * 10.0^(c0 + c1 * logλ + c2 * (logλ * logλ) + c3 * (logλ * logλ * logλ))

        Pe = Korg.kboltz_cgs * (5040.0 / θ) * nₑ
        f * Pe
    end

    @testset "He⁻ free-free absorption" begin
        @testset "compare to Gray polynomial" begin
            T_vals = 2520:100:10_080
            λ_vals = 5063:1000:150_000 # Å
            θ_vals = 5040.0 ./ T_vals
            ν_vals = Korg.c_cgs ./ (λ_vals * 1e-8)

            #assume nₑ = n(He I) = 1
            α_ref = _Heminus_ff_cross_section_Gray.(λ_vals, θ_vals', 1.0)
            #use the unwrapped (not bounds-checked) version for easier broadcasting
            α_korg = Korg.ContinuumAbsorption._Heminus_ff.(ν_vals, T_vals', 1.0, 1.0)

            @test assert_allclose_grid(α_korg, α_ref, [("λ", λ_vals, "Å"), ("T", T_vals, "K")];
                                       rtol=0.05)
        end
        # this really only amounts to a sanity check because the absolute tolerance is of the same
        # magnitude as the actual values
        @testset "Gray (2005) Fig 8.5$panel comparison" for panel in ["b", "c"]
            calculated, ref = Gray_opac_compare.Gray05_comparison_vals(panel, "Heminus_ff")
            @test all(calculated .≥ 0.0)
            @test assert_allclose(calculated, ref; rtol=0, atol=0.02)
        end
    end

    function check_H2plus_ff_and_bf_absorption(target_precision, verbose=true)
        # Table II from Bates (1952) gives the "absorption coefficient" (corrected for stimulated
        # emission) in units of 1e-39 cm⁻¹ / (H atom/cm³) / (H⁺ ion/cm³).
        #
        # Note: we clipped some of the upper rows.  The commented-out columns are temperatures which
        # are currently not supported, but which we would like to support later.
        _table = [1.77 1.35 1.09 0.92 0.79 0.62 0.51;  #  4000
                  1.96 1.47 1.17 0.98 0.84 0.66 0.54;  #  5000
                  2.15 1.57 1.25 1.04 0.89 0.69 0.57;  #  6000
                  2.33 1.67 1.31 1.08 0.92 0.71 0.58;  #  7000
                  2.53 1.77 1.37 1.12 0.95 0.73 0.59;  #  8000
                  2.74 1.87 1.43 1.16 0.97 0.74 0.60;  #  9000
                  2.95 1.97 1.48 1.19 0.99 0.75 0.61;  # 10000
                  3.4 2.18 1.58 1.25 1.03 0.77 0.62;  # 12000
                  4.0 2.40 1.69 1.30 1.06 0.78 0.62;  # 14000
                  4.7 2.64 1.80 1.36 1.09 0.78 0.62;  # 16000
                  5.4 2.91 1.91 1.41 1.11 0.79 0.61;  # 18000
                  6.3 3.2 2.03 1.46 1.14 0.79 0.61;  # 20000
                  7.3 3.5 2.16 1.52 1.16 0.79 0.60;  # 22000
                  8.4 3.8 2.29 1.57 1.18 0.79 0.59;  # 24000
                  9.6 4.2 2.42 1.63 1.21 0.79 0.58]  # 26000
        # the columns of this table are Temperatures (in K)
        _T_vals = [4.0e3, 5.0e3, 6.0e3, 7.0e3, 8.0e3, 1.0e4, 1.2e4] #=2.5e3, 3.0e3, 3.5e3, =#
        # the rows are different wavenumbers (in cm⁻¹)
        _wavenumbers = [4000, 5000, 6000, 7000, 8000, 9000, 10000, 12000, 14000, 16000, 18000,
            20000, 22000, 24000, 26000]

        # to match the implicit assumption in Bates (1952) that ndens(H I) = ndens(H I, n = 1)
        partition_val = 2.0

        # pick stupid dummy values (we're just going to remove them later):
        nH_I = 1e15
        nH_II = 1e13

        const_factor = 1e39 / (nH_I * nH_II)
        function calc_func(ν, T)
            const_factor * Korg.ContinuumAbsorption.H2plus_bf_and_ff([ν], T, nH_I, nH_II)[1]
        end

        λ_cgs = 1.0 ./ _wavenumbers
        _compare_against_table(_table, λ_cgs, _T_vals, calc_func, target_precision, verbose)
    end

    @testset "combined H₂⁺ ff and bf absorption" begin
        @test check_H2plus_ff_and_bf_absorption(0.2)
        @testset "Gray (2005) Fig 8.5$panel comparison" for panel in ["b"]
            calculated, ref = Gray_opac_compare.Gray05_comparison_vals(panel, "H2plus")
            @test all(calculated .≥ 0.0)
            @test assert_allclose(calculated, ref; rtol=0, atol=0.02)
        end
    end

    """
        compare_gauntff_kurucz([rtol])

    Compare the free-free gaunt factor to the reference values tabulated in section 5.1 of Kurucz
    (1970). These reference tabulated values are less accurate than the values we're actually using
    (since these values were derived from a figure in Karsas and Latter (1961)), but they're adequate
    for testing purposes.
    """
    function compare_gauntff_kurucz(rtol=0.15)

        # First, get bounds of log₁₀(u) = log₁₀(RydbergH*Z²/(k*T)) & log₁₀(γ²) = log₁₀(Rydberg*Z²/(k*T)),
        # in which our actual data is known.
        _log10_u_bounds, _log10_γ2_bounds = bounds(Korg.ContinuumAbsorption._gauntff_interpolator.itp)
        log10_u_bounds = Korg.Interval(_log10_u_bounds...)
        log10_γ2_bounds = Korg.Interval(_log10_γ2_bounds...)

        # Next, initialize the reference data
        log10_γ2 = [-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]
        log10_u = [-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5]
        ref_gaunt_ff = [5.53 5.49 5.46 5.43 5.40 5.25 5.00 4.69 4.48 4.16 3.85;
                        4.91 4.87 4.84 4.80 4.77 4.63 4.40 4.13 3.87 3.52 3.27;
                        4.29 4.25 4.22 4.18 4.15 4.02 3.80 3.57 3.27 2.98 2.70;
                        3.64 3.61 3.59 3.56 3.54 3.41 3.22 2.97 2.70 2.45 2.20;
                        3.00 2.98 2.97 2.95 2.94 2.81 2.65 2.44 2.21 2.01 1.81;
                        2.41 2.41 2.41 2.41 2.41 2.32 2.19 2.02 1.84 1.67 1.50;
                        1.87 1.89 1.91 1.93 1.95 1.90 1.80 1.68 1.52 1.41 1.30;
                        1.33 1.39 1.44 1.49 1.55 1.56 1.51 1.42 1.33 1.25 1.17;
                        0.90 0.95 1.00 1.08 1.17 1.30 1.32 1.30 1.20 1.15 1.11;
                        0.45 0.48 0.52 0.60 0.75 0.91 1.15 1.18 1.15 1.11 1.08;
                        0.33 0.36 0.39 0.46 0.59 0.76 0.97 1.09 1.13 1.10 1.08;
                        0.19 0.21 0.24 0.28 0.38 0.53 0.76 0.96 1.08 1.09 1.09]

        # now, get the set of log10_γ2 and log10_u values where we can make the comparison
        ref_u_slc = Korg.contained_slice(log10_u, log10_u_bounds)
        ref_γ2_slc = Korg.contained_slice(log10_γ2, log10_γ2_bounds)
        @assert length(ref_u_slc) > 0 && length(ref_γ2_slc) > 0
        cmp_log10_u, cmp_log10_γ2 = log10_u[ref_u_slc], log10_γ2[ref_γ2_slc]

        # finally, perform the comparison
        assert_allclose_grid(Korg.ContinuumAbsorption._gauntff_interpolator.(cmp_log10_u,
                                                                             cmp_log10_γ2'),
                             ref_gaunt_ff[ref_u_slc, ref_γ2_slc],
                             [("log₁₀u", cmp_log10_u), ("log₁₀γ²", cmp_log10_γ2)];
                             rtol=rtol, atol=0.0)
    end

    """
        gray_H_I_ff_absorption_coef(λ, T, [gaunt_func])

    Computes the H I free-free linear absorption coefficient per neutral Hydrogen atom by adapting
    equation 8.10 from Gray (2005). This includes the correction for stimulated emission under LTE.
    In contrast to `Korg.ContinuumAbsorption.positive_ion_ff`, this implementation directly combines the Saha
    equation and uses constants that have been cast into a completely separate format.

    # Arguments

      - `λ`: wavelength in Å
      - `T`: temperature in K
      - `gaunt_func`: optional argument specifying a function to compute the free-free gaunt factor. The
        function should accept log_u and log_γ2 as arguments. By default, this uses the approximation
        provided by equation 8.6 from Gray (2005).

    # Notes

    I worked out the algebra in detail and can confirm that the this is assumes that partition function
    of H I is 2.0 at all temperatures.
    """
    function gray_H_I_ff_absorption_coef(λ, T, gaunt_func=nothing)

        # adapt equation 8.10
        R = 1.0968e5 # cm⁻¹
        I = Korg.hplanck_eV * Korg.c_cgs * R # eV
        θ = 5040 / T # eV^-1
        loge = log10(MathConstants.e)
        α₀ = 1.0449e-26 #cm² Å⁻³
        χ_λ = 1.2398e4 / λ # eV
        Z2 = 1 # this is an implicit assumption from the textbook

        gaunt_ff = if isnothing(gaunt_func)
            # equation 8.6 from Gray (2005)
            temp_R = R / 1e8 # Å⁻¹
            1.0 + 0.3456 * (0.5 + loge / (θ * χ_λ)) / cbrt(λ * temp_R)
        else
            β = 1.0 / (Korg.kboltz_eV * T)
            ν = Korg.c_cgs * 1e8 / λ
            log_u = log10(Korg.hplanck_eV * ν * β)
            log_γ2 = log10(Korg.RydbergH_eV * Z2 * β)
            gaunt_func(log_u, log_γ2)
        end

        uncorrected = α₀ * λ^3 * gaunt_ff * loge * 10.0^(-θ * I) / (2 * θ * I)
        # according to equation 8.18, we need to correct the value for stimulated emission
        uncorrected * (1 - 10.0^(-χ_λ * θ))
    end

    @testset "H I free-free absorption" begin
        @testset "Kurucz (1970) Free-Free Gaunt Factor comparison" begin
            @test compare_gauntff_kurucz(0.15)
        end

        @testset "Gray (2005) H I ff implementation comparison" begin
            # Compare the different formulations of the free-free absorption under the conditions of
            # the panels of figure 8.5 of Gray (2005). Outside of these conditions, it's unclear where
            # Gray's approximation for the gaunt factor breaks down.
            # This comparison is much more constraining than directly comparing against the relevant
            # curve extracted from the panel because the relevant curve includes both free-free and
            # bound-free absorption. It's also much more constraining when H I ff absorption is
            # subdominant
            χs = [(Korg.RydbergH_eV, -1.0, -1.0)]
            # implicitly assumed by the Gray implementation
            Us = Dict([Korg.species"H_I" => (T -> 2.0), Korg.species"H_II" => (T -> 1.0)])

            λ_vals = [3e3 + (i - 1) * 250 for i in 1:69] # equally spaced vals from 3e3 Å through 2e4 Å
            ν_vals = Korg.c_cgs * 1e8 ./ λ_vals

            for (logPₑ, T) in [(1.08, 5143.0), (1.77, 6429.0), (2.50, 7715.0), (2.76, 11572.0)]
                ref_absorption_coef = gray_H_I_ff_absorption_coef.(λ_vals, T)

                nₑ = 10.0^logPₑ / (Korg.kboltz_cgs * T)
                nH_total = nₑ * 100.0 # this is totally arbitrary

                wII, wIII = Korg.saha_ion_weights(T, nₑ, 1, χs, Us)
                nH_I = nH_total / (1 + wII)
                nH_II = nH_total * wII / (1 + wII)
                # recall that H I ff actually refers to: photon + e⁻ + H II -> e⁻ + H II
                absorption_coef = zeros(length(ν_vals))
                Korg.ContinuumAbsorption.positive_ion_ff_absorption!(absorption_coef, ν_vals, T,
                                                                     Dict([Korg.species"H II" => nH_II]),
                                                                     nₑ)
                @test assert_allclose_grid(absorption_coef ./ nH_I, ref_absorption_coef,
                                           [("λ", λ_vals, "Å")];
                                           rtol=0.034, atol=0)
                @test all(absorption_coef .≥ 0.0)
            end
        end
    end

    """
        gray_H_I_bf_absorption_coef(λ, T)

    Computes the H I bound-free linear absorption coefficient per neutral Hydrogen atom (with the
    stimulated emission correction term) by adapting equation 8.8 from Gray (2005).

    In contrast to `Korg.ContinuumAbsorption.H_I_bf`, this implementation:

      - uses a completely different approximation for computing the cross-section (this includes
        differences in gaunt factor estimation).
      - uses different logic for approximating the contributions of higher energy levels to the
        absorption. In this approximation the lowest 2 energy levels that can be ionized at a given
        wavelength, λ, are always explicitly summed and the contributions from all higher energy levels
        are always integrated.
      - uses constants that have been cast into a completely different format.

    # Arguments

      - `λ`: wavelength in Å
      - `T`: temperature in K

    # Notes

    This implicitly assumes that the H I partition function is always equal to 2.    # equation 8.8 gives the formula without stimulated emission corrections
    """
    function gray_H_I_bf_absorption_coef(λ, T)
        # equation 8.8 gives the formula without stimulated emission corrections
        # λ must have units of Å and T must have units of K

        θ = 5040 / T # eV^-1
        photon_energy_eV = Korg.hplanck_eV * Korg.c_cgs * 1e8 / λ

        # They are not particularly tranparent about this, but I believe n0 should be the lowest energy
        # state for which the ionization energy is greater than or equal to photon_energy_eV
        full_ion_energy = Korg.ionization_energies[1][1] # ionization from ground state

        # to ionize electron in state n, need:   photon_energy_eV ≥ full_ion_energy /n²
        # Rearrange eqn to determine values of n that can be excited by a given photon_energy_eV:
        #                                                      n  ≥ √(full_ion_energy/photon_energy_eV)
        # Then n₀ is the minimum value of n that satisfies this:
        n₀ = ceil(Int64, sqrt(full_ion_energy / photon_energy_eV))

        # Let's compute the summation in equation 8.8 that runs from n₀ through n₀+2
        Σ_term = 0.0
        for n in n₀:(n₀+3)
            # compute excitation potential, χ, (from equation 8.3):
            χ = full_ion_energy * (1.0 - 1.0 / (n * n)) #eV
            # estimate bound-free gaunt factor (from equation 8.5):
            R = 1.0968e-3 # Å⁻¹
            gaunt_bf = 1.0 - 0.3456 * (λ * R / (n * n) - 0.5) / cbrt(λ * R)

            Σ_term += 10.0^(-θ * χ) * gaunt_bf / n^3
        end

        # Let's compute the "integal_term"
        χ₃ = full_ion_energy * (1 - 1 / (n₀ + 3)^2)
        loge = 0.43429 # = log₁₀e
        integral_term = (10^(-χ₃ * θ) - 10^(-full_ion_energy * θ)) * 0.5 * loge /
                        (θ * full_ion_energy)

        α₀ = 1.0449e-26 #cm² Å⁻³
        uncorrected_absorption_coef = α₀ * λ^3 * (Σ_term + integral_term)

        # according to equation 8.18, we need to correct the value for stimulated emission
        χ_λ = 1.2398e4 / λ
        uncorrected_absorption_coef * (1 - 10.0^(-χ_λ * θ))
    end

    @testset "H I bound-free absorption" begin
        @testset "Gray (2005) implementation comparison" begin
            # Compare the different formulations of the bound-free absorption under the conditions of
            # the panels of figure 8.5 of Gray (2005). Outside of these conditions, it's unclear
            # if/where Gray's approximation for the gaunt factor breaks down.
            # This comparison is much more constraining than directly comparing against the relevant
            # curve extracted from the panel because the relevant curve includes both free-free and
            # bound-free absorption. It's also much more constraining when H I bf absorption is
            # subdominant

            H_I_partition_val = 2.0 # implicitly assumed by the Gray implementation

            λ_vals = [3e3 + (i - 1) * 250 for i in 1:69] # equally spaced vals from 3e3 Å through 2e4 Å
            ν_vals = Korg.c_cgs * 1e8 ./ λ_vals

            for T in [5143.0, 6429.0, 7715.0, 11572.0]
                ref_absorption_coef = gray_H_I_bf_absorption_coef.(λ_vals, T)

                absorption_coef = reverse(Korg.ContinuumAbsorption.H_I_bf(reverse(ν_vals), T, 1, 0,
                                                                          0,
                                                                          1 / H_I_partition_val))
                @test assert_allclose_grid(absorption_coef, ref_absorption_coef,
                                           [("λ", λ_vals, "Å")];
                                           rtol=0.09, print_rachet_info=false)
                @test all(absorption_coef .≥ 0.0)
            end
        end
    end

    # this is only included here for completeness.
    # Given the strong consistency when comparing the bf and ff values individually against the
    # alternative implementations described in Gray (2005), I'm unconcerned by the need for large
    # tolerances.
    @testset "combined H I bf and ff absorption" begin
        @testset "Gray (2005) Fig 8.5$panel comparison" for (panel, atol) in [("b", 0.04),
            ("c", 0.125),
            ("d", 30.0)]
            calculated, ref = Gray_opac_compare.Gray05_comparison_vals(panel, "H")
            @test all(calculated .≥ 0.0)
            @test assert_allclose(calculated, ref; rtol=0, atol=atol)
        end
    end

    @testset "TOPbase/NORAD bound-free absorption" begin
        T = 7800.0 #K, this is fairly arbitrary
        ndens_species = 3.0e16 #cm⁻³, this is fairly arbitrary

        λ_vals = OP_compare._dflt_λ_vals(Korg.species"H I")

        # compute the absorption coefficients using the TOPbase/NORAD data
        H_I_bf_OP = OP_compare.calc_hydrogenic_bf_absorption_coef(λ_vals, T, ndens_species,
                                                                  Korg.species"H I";
                                                                  use_OP_data=true)
        # compute the absorption coefficients using our function for hydrogenic atoms
        H_I_bf_default = OP_compare.calc_hydrogenic_bf_absorption_coef(λ_vals, T, ndens_species,
                                                                       Korg.species"H I";
                                                                       use_OP_data=false)

        λ_comp_intervals = [(λ_vals[1], λ_vals[end])]
        @test assert_allclose_grid(H_I_bf_OP, H_I_bf_default, [("λ", λ_vals, "Å")];
                                   rtol=OP_compare._hydrogenic_rtol(Korg.species"H I"), atol=0.0,
                                   err_msg=("\nH I bf absorption coefficients " *
                                            "computed using data from the opacity project are " *
                                            "inconsistent with the default implementation."),
                                   print_rachet_info=true)
    end
end
