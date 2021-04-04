include("opacity_comparison_funcs.jl")

function _calc_Hminus_ff_absorption_coef(ν, T)
    # We invert the H⁻ opacity to solve for the absorption coefficient:
    # α_ff(H⁻) = κ_ν * ρ / (n(H I) * Pₑ)
    # See the description of Hminus_ff for more details

    # The choices of nH_I, ne, and ρ are completely unimportant since α_ff(H⁻) has no explicit
    # dependence on any of them.
    nH_I, ne = 1e17, 1e14
    ρ = 1e-24 * nH_I # the value is a little nonsensical but it's the right order of magnitude

    # To allow the best comparison to real data lets set the H I partition function to 2, so that
    # n(H I, n = 1) = n(H I).
    nH_I_div_partition = nH_I/2.0

    κ_ν = SSSynth.ContinuumOpacity.Hminus_ff(nH_I_div_partition, ne, ν, ρ, T)
    Pₑ = ne * SSSynth.kboltz_cgs * T
    αff_H⁻ = κ_ν * ρ / (nH_I * Pₑ)

    # finally divide by 1e-26 to match the units of the table
    αff_H⁻/1e-26
end


function check_Hminus_ff_values(target_precision = 0.01, verbose = true)
    # I tabulated this data from Bell & Berrington (1987)
    # each row corresponds to a separate λ (in Å)
    _ff_λ_vals = [151890, 113918, 91134, 45567, 30378, 22784, 18227, 15189, 13019, 11392, 10126,
                  9113, 7595, 6510, 5696, 5063, 4557, 3645, 3038, 2604, 2278, 1823]
    # each column corresponds to a separate θ value where θ = 5040/T
    _ff_θ_vals = [0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.8, 3.6]
    _ff_T_vals = 5040.0./_ff_θ_vals
    _ff_table = [75.1   90     118    144    168    191    212    234    253    325    388;
                 42.3   50.6   66.4   80.8   94.5   107    120    131    142    183    219;
                 27     32.4   42.6   51.9   60.7   68.9   76.8   84.2   91.4   117    140;
                 6.79   8.16   10.7   13.1   15.3   17.4   19.4   21.2   23     29.5   35;
                 3.04   3.65   4.8    5.86   6.86   7.79   8.67   9.5    10.3   13.2   15.6;
                 1.73   2.08   2.74   3.37   3.9    4.5    5.01   5.5    5.95   7.59   9.06;
                 1.11   1.34   1.77   2.17   2.53   2.87   3.2    3.51   3.8    4.92   5.97;
                 0.781  0.94   1.24   1.52   1.78   2.02   2.26   2.48   2.69   3.52   4.31;
                 0.579  0.699  0.924  1.13   1.33   1.51   1.69   1.86   2.02   2.67   3.31;
                 0.448  0.539  0.711  0.871  1.02   1.16   1.29   1.43   1.57   2.09   2.6;
                 0.358  0.432  0.572  0.702  0.825  0.943  1.06   1.17   1.28   1.73   2.17;
                 0.293  0.354  0.468  0.576  0.677  0.777  0.874  0.969  1.06   1.45   1.83;
                 0.208  0.25   0.332  0.409  0.484  0.557  0.63   0.702  0.774  1.06   1.36;
                 0.154  0.188  0.249  0.309  0.367  0.424  0.482  0.539  0.597  0.83   1.07;
                 0.121  0.146  0.195  0.241  0.288  0.334  0.381  0.428  0.475  0.667  0.861;
                 0.0965 0.117  0.157  0.195  0.234  0.272  0.311  0.351  0.39   0.549  0.711;
                 0.0791 0.0959 0.129  0.161  0.194  0.227  0.26   0.293  0.327  0.463  0.602;
                 0.052  0.0633 0.0859 0.108  0.131  0.154  0.178  0.201  0.225  0.321  0.418;
                 0.0364 0.0447 0.0616 0.0789 0.0966 0.114  0.132  0.15   0.169  0.243  0.318;
                 0.0277 0.0342 0.0476 0.0615 0.076  0.0908 0.105  0.121  0.136  0.199  0.262;
                 0.0228 0.028  0.0388 0.0499 0.0614 0.0732 0.0851 0.0972 0.11   0.16   0.211;
                 0.0178 0.0222 0.0308 0.0402 0.0498 0.0596 0.0695 0.0795 0.0896 0.131  0.172]

    # only look at the subtable that might be relevant
    sub_table = view(_ff_table, 2:20, 1:9)
    comp_λ = view(_ff_λ_vals, 2:20)
    comp_T = view(_ff_T_vals, 1:9)

    success = true
    ν_vals = (SSSynth.c_cgs*1e8)./comp_λ
    for index = 1:length(comp_T)
        cur_T = comp_T[index]
        ref_coef = view(sub_table, :, index)
        calc_coef = _calc_Hminus_ff_absorption_coef.(ν_vals, cur_T)
        precision = (abs.(calc_coef .- ref_coef)./ref_coef)
        max_err, max_err_ind = findmax(precision)
        if max_err > target_precision
            if verbose
                if success
                    println("There is a problem. The max error of polynomial should be: ",
                            target_precision*100, "%.")
                end
                println("Max error for T = ", cur_T, "K is: ", 100 * max_err,
                        "% at λ = ", comp_λ[max_err_ind], " Å.")
            end
            success = false
        end
    end
    return success
end

@testset "H⁻ free-free opacity" begin
    @test check_Hminus_ff_values(0.0225)
    # we only have measurements for b and c
    @testset "Gray (2005) Fig 8.5$panel comparison" for panel in ["b", "c"]
        calculated, ref = Gray05_comparison_vals(panel,"Hminus_ff")
        @test all(calculated .≥ 0.0)
        @test all(abs.(calculated - ref) .≤ Gray05_atols[panel])
    end
end

function check_Hminus_bf_values(target_precision = 0.002, verbose = true)
    # this is the tabulated data from Wishart (1979)
    _bf_table = [16300 0.1989; 16200 0.4974; 16100 0.8697; 16000 1.302; 15750 2.575; 15500 4.052;
                 15250 5.677;  15000 7.407;  14750 9.211;  14500 11.07; 14250 12.95; 14000 14.85;
                 13750 16.74;  13500 18.62;  13250 20.46;  13000 22.26; 12750 24.02; 12500 25.71;
                 12250 27.33;  12000 28.87;  11750 30.34;  11500 31.72; 11250 33.01; 11000 34.19;
                 10750 35.28;  10500 36.25;  10250 37.13;  10000 37.89; 9750 38.53;  9500 39.06;
                 9250 39.48;   9000 39.77;   8750 39.95;   8500 40.01;  8250 39.95;  8000 39.77;
                 7750 39.48;   7500 39.07;   7250 38.54;   7000 37.91;  6750 37.17;  6500 36.32;
                 6250 35.37;   6000 34.32;   5750 33.17;   5500 31.94;  5250 30.62;  5000 29.23;
                 4750 27.77;   4500 26.24;   4250 24.65;   4000 23.02;  3750 21.35;  3500 19.65;
                 3250 17.92;   3000 16.19;   2750 14.46;   2500 12.75;  2250 11.08;  2000 9.453;
                 1750 7.918;   1500 6.512;   1250 5.431]

    _bf_λ_vals = view(_bf_table, :, 1) # in Ångstroms
    _bf_α_vals = view(_bf_table, :, 2) .* 1e-18 # in cm² per H⁻ particle.
    comp_λ_vals = view(_bf_λ_vals, 8:59)
    comp_α_vals = view(_bf_α_vals, 8:59)
    ν_vals = (SSSynth.c_cgs*1e8)./comp_λ_vals

    coefs = SSSynth.ContinuumOpacity._Hminus_bf_cross_section.(ν_vals)
    precision = abs.(coefs .- comp_α_vals)./comp_α_vals
    max_err, max_err_ind = findmax(precision)

    if max_err < target_precision
        return true
    else
        if verbose
            println("There is a problem. The max error of polynomial should be: ",
                    target_precision*100, "%.")
            println("Max error is: ", 100 * max_err,
                    "% at λ = ", comp_λ_vals[max_err_ind], " Å.")
        end
        return false
    end
end


@testset "H⁻ bound-free opacity" begin
    @test check_Hminus_bf_values(0.0025)
    @testset "Gray (2005) Fig 8.5$panel comparison" for panel in ["a", "b", "c"]
        calculated, ref = Gray05_comparison_vals(panel,"Hminus_bf")
        @test all(calculated .≥ 0.0)
        @test all(abs.(calculated - ref) .≤ Gray05_atols[panel])
    end
end


@testset "He⁻ free-free opacity" begin
    # this really only amounts to a sanity check because the absolute tolerance is of the same
    # magnitude as the actual values
    @testset "Gray (2005) Fig 8.5$panel comparison" for panel in ["b", "c"]
        calculated, ref = Gray05_comparison_vals(panel,"Heminus_ff")
        @test all(calculated .≥ 0.0)
        @test all(abs.(calculated - ref) .≤ Gray05_atols[panel])
    end
end

function check_H2plus_ff_and_bf_opacity(target_precision, verbose = true)
    # table II from Bates (1952) gives the absorption coefficient (corrected for stimulated
    # emission) in units of 1e-39 cm⁻¹ / (H atom/cm³) / (H⁺ ion/cm³).
    # Note: we clipped some of the upper rows
    _table = [  3.4  2.60  2.10  1.77  1.35  1.09  0.92  0.79  0.62  0.51;  #  4000
                4.0  2.98  2.36  1.96  1.47  1.17  0.98  0.84  0.66  0.54;  #  5000
                4.8  3.4   2.63  2.15  1.57  1.25  1.04  0.89  0.69  0.57;  #  6000
                5.6  3.9   2.91  2.33  1.67  1.31  1.08  0.92  0.71  0.58;  #  7000
                6.7  4.4   3.2   2.53  1.77  1.37  1.12  0.95  0.73  0.59;  #  8000
                7.9  5.0   3.5   2.74  1.87  1.43  1.16  0.97  0.74  0.60;  #  9000
                9.3  5.6   3.9   2.95  1.97  1.48  1.19  0.99  0.75  0.61;  # 10000
               13.0  7.2   4.7   3.4   2.18  1.58  1.25  1.03  0.77  0.62;  # 12000
               18.1  9.3   5.7   4.0   2.40  1.69  1.30  1.06  0.78  0.62;  # 14000
               25.1 11.9   7.0   4.7   2.64  1.80  1.36  1.09  0.78  0.62;  # 16000
               35   15.2   8.4   5.4   2.91  1.91  1.41  1.11  0.79  0.61;  # 18000
               47   19.3  10.2   6.3   3.2   2.03  1.46  1.14  0.79  0.61;  # 20000
               64   24.3  12.2   7.3   3.5   2.16  1.52  1.16  0.79  0.60;  # 22000
               86   31    14.6   8.4   3.8   2.29  1.57  1.18  0.79  0.59;  # 24000
              114   38    17.3   9.6   4.2   2.42  1.63  1.21  0.79  0.58]  # 26000
    # the columns of this table are Temperatures (in K)
    _T_vals = [2.5e3 3.0e3 3.5e3 4.0e3 5.0e3 6.0e3 7.0e3 8.0e3 1.0e4 1.2e4]
    # the rows are different wavenumbers (in cm⁻¹)
    _wavenumbers = [ 4000,  5000,  6000,  7000,  8000,  9000, 10000, 12000, 14000, 16000, 18000,
                    20000, 22000, 24000, 26000]

    sub_table = _table
    comp_λ = 1.0e8 ./ _wavenumbers
    comp_T = _T_vals
    ν_vals = _wavenumbers .* SSSynth.c_cgs

    # to match the implicit assumption in Bates (1952) that ndens(H I) = ndens(H I, n = 1)
    partition_val = 2.0

    # pick stupid dummy values (we're just going to remove them later):
    nH_I = 1e15
    nH_II = 1e13
    ρ = 1.0

    nH_I_divU = nH_I/partition_val
    success = true
    for index = 1:length(comp_T)
        cur_T = comp_T[index]
        ref_coef = view(sub_table, :, index)
        calc_opacity = SSSynth.ContinuumOpacity.H2plus_bf_and_ff.(nH_I_divU, nH_II, ν_vals, ρ,
                                                                  cur_T)
        calc_coef = 1e39 * calc_opacity * ρ/(nH_I * nH_II)
        precision = (abs.(calc_coef .- ref_coef)./ref_coef)
        max_err, max_err_ind = findmax(precision)
        if max_err > target_precision
            if verbose
                if success
                    println("There is a problem. The max error of polynomial should be: ",
                            target_precision*100, "%.")
                end
                println("Max error for T = ", cur_T, "K is: ", 100 * max_err,
                        "% at λ = ", comp_λ[max_err_ind], " Å.")
            end
            success = false
        end
    end
    return success
end

@testset "combined H₂⁺ ff and bf opacity" begin
    @test check_H2plus_ff_and_bf_opacity(0.015)
    @testset "Gray (2005) Fig 8.5$panel comparison" for panel in ["b"]
        calculated, ref = Gray05_comparison_vals(panel,"H2plus")
        @test all(calculated .≥ 0.0)
        @test all(abs.(calculated - ref) .≤ Gray05_atols[panel])
    end
end


"""
    gray_H_I_ff_absorption_coef(λ, T, [gaunt_func])

Computes the H I free-free absorption coefficient per neutral Hydrogen atom by adapting equation
8.10 from Gray (2005). In contrast to `SSSynth.ContinuumOpacity.H_I_ff`, this implementation
directly combines the Saha equation and uses constants that have been cast into a completely
separate format.

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
function gray_H_I_ff_absorption_coef(λ, T, gaunt_func = nothing)

    # adapt equation 8.10
    R = 1.0968e5 # cm⁻¹
    I = SSSynth.hplanck_eV * SSSynth.c_cgs * R # eV
    θ = 5040/T # eV^-1
    loge = log10(MathConstants.e)
    α₀ = 1.0449e-26 #cm² Å⁻³
    χ_λ = 1.2398e4/λ # eV
    Z2 = 1 # this is an implicit assumption from the textbook

    gaunt_ff = if isnothing(gaunt_func)
        # equation 8.6 from Gray (2005)
        temp_R = R/1e8 # Å⁻¹
        1.0 + 0.3456*(0.5 + loge/(θ*χ_λ))/cbrt(λ*temp_R)
    else
        β = 1.0/(SSSynth.kboltz_eV * T)
        ν = SSSynth.c_cgs * 1e8 / λ
        log_u = log10(SSSynth.hplanck_eV * ν * β)
        log_γ2 = log10(SSSynth.RydbergH_eV * Z2 * β)
        gaunt_func(log_u, log_γ2)
    end

    uncorrected = α₀ *λ^3 * gaunt_ff * loge * 10.0^(-θ*I)/ (2*θ*I)
    # according to equation 8.18, we need to correct the value for stimulated emission
    uncorrected*(1 - 10.0^(-χ_λ * θ))
end


@testset "H I free-free opacity" begin
    @testset "Gray (2005) implementation comparison" begin
        # Compare the different formulations of the free-free opacity under the conditions of the
        # panels of figure 8.5 of Gray (2005). Outside of these conditions, it's unclear where
        # Gray's approximation for the gaunt factor breaks down.
        # This comparison is much more constraining than directly comparing against the relevant
        # curve extracted from the panel because the relevant curve includes both free-free and
        # bound-free absorption. It's also much more constraining when H I ff opacity is subdominant
        χs = [SSSynth.RydbergH_eV, -1.0]
        Us = [T -> 2.0, T -> 1.0] # implicitly assumed by the Gray implementation

        λ_vals = [3e3 + (i-1)*250 for i = 1:69] # equally spaced vals from 3e3 Å through 2e4 Å
        ν_vals = SSSynth.c_cgs*1e8./λ_vals

        for (logPₑ, T) in [(1.08, 5143.0), (1.77, 6429.0), (2.50, 7715.0), (2.76, 11572.0)]
            ref_absorption_coef = gray_H_I_ff_absorption_coef.(λ_vals, T)

            nₑ = 10.0^logPₑ / (SSSynth.kboltz_cgs * T)
            nH_total = nₑ * 100.0 # this is totally arbitrary
            ρ = nH_total * 1.67e-24/0.76 # this is totally arbitrary

            weights = SSSynth.saha(χs, Us, T, nₑ)
            nH_I = nH_total * weights[1]
            nH_II = nH_total * weights[2]
            # recall that H I ff actually refers to: photon + e⁻ + H II -> e⁻ + H II
            absorption_coef = SSSynth.ContinuumOpacity.H_I_ff.(nH_II, nₑ, ν_vals, ρ, T) * ρ/nH_I
            @test maximum(abs.(absorption_coef - ref_absorption_coef)/absorption_coef) < 0.0015
            @test all(absorption_coef .≥ 0.0)
        end
    end
end


"""
    gray_H_I_bf_absorption_coef(λ, T)

Computes the H I bound-free absorption coefficient per neutral Hydrogen atom (with the stimulated
emission correction term) by adapting equation 8.8 from Gray (2005). 

In contrast to `SSSynth.ContinuumOpacity.H_I_bf`, this implementation:
- uses a completely different approximation for computing the cross-section (this includes
  differences in gaunt factor estimation).
- uses different logic for approximating the contributions of higher energy levels to the opacity. 
  In this approximation the lowest 2 energy levels that can be ionized at a given wavelength, λ,
  are always explicitly summed and the contributions from all higher energy levels are always 
  integrated.
- uses constants that have been cast into a completely different format.

# Arguments
- `λ`: wavelength in Å
- `T`: temperature in K

# Notes
This implicitly assumes that the H I partition function is always equal to 2.
"""
function gray_H_I_bf_absorption_coef(λ, T)
    # equation 8.8 gives the formula without stimulated emission corrections
    # λ must have units of Å and T must have units of K

    θ = 5040/T # eV^-1
    photon_energy_eV = SSSynth.hplanck_eV * SSSynth.c_cgs * 1e8/λ

    # They are not particularly tranparent about this, but I believe n0 should be the lowest energy
    # state for which the ionization energy is greater than or equal to photon_energy_eV
    full_ion_energy = 13.6 # ionization from ground state

    # to ionize electron in state n, need:   photon_energy_eV ≥ full_ion_energy /n²
    # Rearrange eqn to determine values of n that can be excited by a given photon_energy_eV:
    #                                                      n  ≥ √(full_ion_energy/photon_energy_eV)
    # Then n₀ is the minimum value of n that satisfies this:
    n₀ = ceil(Int64, sqrt(full_ion_energy / photon_energy_eV))

    # Let's compute the summation in equation 8.8 that runs from n₀ through n₀+2
    Σ_term = 0.0
    for n = n₀:(n₀+2)
        # compute excitation potential, χ, (from equation 8.3):
        χ = full_ion_energy * (1.0 - 1.0/(n*n)) #eV
        # estimate bound-free gaunt factor (from equation 8.5):
        R = 1.0968e-3 # Å⁻¹
        gaunt_bf = 1.0 - 0.3456*(λ*R/(n*n) - 0.5)/cbrt(λ*R)

        Σ_term += 10.0^(-θ*χ) * gaunt_bf / n^3
    end

    # Let's compute the "integal_term"
    χ₃ = full_ion_energy * (1 - 1/(n₀+3)^2)
    loge = 0.43429 # = log₁₀e
    integral_term = (10^(-χ₃*θ) - 10^(-full_ion_energy*θ)) * 0.5 * loge / (θ * full_ion_energy)

    α₀ = 1.0449e-26 #cm² Å⁻³
    uncorrected_absorption_coef = α₀ * λ^3 * (Σ_term + integral_term)

    # according to equation 8.18, we need to correct the value for stimulated emission
    χ_λ = 1.2398e4/λ
    uncorrected_absorption_coef * (1 - 10.0^(-χ_λ * θ))
end

@testset "H I bound-free opacity" begin
    @testset "Gray (2005) implementation comparison" begin
        # Compare the different formulations of the bound-free opacity under the conditions of the
        # panels of figure 8.5 of Gray (2005). Outside of these conditions, it's unclear if/where
        # Gray's approximation for the gaunt factor breaks down.
        # This comparison is much more constraining than directly comparing against the relevant
        # curve extracted from the panel because the relevant curve includes both free-free and
        # bound-free absorption. It's also much more constraining when H I bf opacity is subdominant

        H_I_ion_energy = SSSynth.RydbergH_eV # use this to decouple test from ionization energies
        H_I_partition_val = 2.0 # implicitly assumed by the Gray implementation

        λ_vals = [3e3 + (i-1)*250 for i = 1:69] # equally spaced vals from 3e3 Å through 2e4 Å
        ν_vals = SSSynth.c_cgs*1e8./λ_vals

        for (logPₑ, T) in [(1.08, 5143.0), (1.77, 6429.0), (2.50, 7715.0), (2.76, 11572.0)]
            ref_absorption_coef = gray_H_I_bf_absorption_coef.(λ_vals, T)

            nₑ = 10.0^logPₑ / (SSSynth.kboltz_cgs * T)
            nH_I = nₑ * 100.0 # this is totally arbitrary
            ρ = nH_I * 1.67e-24/0.76 # this is totally arbitrary

            absorption_coef = SSSynth.ContinuumOpacity.H_I_bf.(nH_I/H_I_partition_val, ν_vals, ρ,
                                                               T, H_I_ion_energy) * ρ/nH_I
            #println(maximum(abs.(absorption_coef - ref_absorption_coef)/absorption_coef))
            @test maximum(abs.(absorption_coef - ref_absorption_coef)/absorption_coef) < 0.004
            @test all(absorption_coef .≥ 0.0)
        end
    end

    @testset "Integral-Summation Equivalence" begin
        # SSSynth.ContinuumOpacity.H_I_bf approximates the sum of the absorption contributions from
        # bound-free transitions originating from high energy levels with an integral. These tests
        # check that the integral does a reasonable job reproducing the direct sum.

        H_I_ion_energy = SSSynth.RydbergH_eV # use this to decouple test from ionization energies
        H_I_partition_val = 2.0 # implicitly assumed by the Gray implementation

        λ_vals = [3e3 + (i-1)*500 for i = 1:35] # equally spaced vals from 3e3 Å through 2e4 Å
        ν_vals = SSSynth.c_cgs*1e8./λ_vals

        nstop_sum = 100

        for (logPₑ, T) in [(1.08, 5143.0), (1.77, 6429.0), (2.50, 7715.0), (2.76, 11572.0)]
            ref_absorption_coef = gray_H_I_bf_absorption_coef.(λ_vals, T)

            nₑ = 10.0^logPₑ / (SSSynth.kboltz_cgs * T)
            nH_I = nₑ * 100.0 # this is totally arbitrary
            ρ = nH_I * 1.67e-24/0.76 # this is totally arbitrary

            opacity_integral = SSSynth.ContinuumOpacity.H_I_bf.(nH_I/H_I_partition_val, ν_vals, ρ,
                                                                T, H_I_ion_energy)
            opacity_sum = SSSynth.ContinuumOpacity.H_I_bf.(nH_I/H_I_partition_val, ν_vals, ρ,
                                                           T, H_I_ion_energy, nstop_sum, false)
            #println(maximum(abs.(opacity_integral - opacity_sum)/opacity_sum))
            @test maximum(abs.(opacity_integral - opacity_sum)/opacity_sum) < 0.003
        end
    end
end


# this is only included here for completeness.
# Given the strong consistency when comparing the bf and ff values individually against the
# alternative implementations described in Gray (2005), I'm unconcerned by the need for large
# tolerances.
@testset "combined H I bf and ff opacity" begin
    @testset "Gray (2005) Fig 8.5$panel comparison" for (panel, atol) in [("b",0.035),
                                                                          ("c",0.125),
                                                                          ("d",35)]
        calculated, ref = Gray05_comparison_vals(panel,"H")
        @test all(calculated .≥ 0.0)
        @test all(abs.(calculated - ref) .≤ atol)
    end
end
