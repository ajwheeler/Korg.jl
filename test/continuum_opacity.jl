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
