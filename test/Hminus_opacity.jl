"""
    calc_H⁻_ff_absorption_coef()

Computes the atomic free-free absorption coefficient, α for H⁻ in the form of equation (8.13) from 
Gray (2005). This will have units of 1e-26 cm^2 per H atom per unit Pₑ (or 1e-26 cm^2 per H atom 
per dynes/cm²). The absorption coefficient includes the correction from stimulated emission.
Note that Pₑ is the partial pressure contribution from the free electrons or Pₑ = nₑ*k*T.

This is function primarily exists to facillitate comparisons against the tabulated data from Bell &
Berrington (1987), who provide their data in this same form.

# Notes
From equations 8.13, 8.18, and 8.19 the free-free opacity from H⁻ is given by:
                              n(H I)          n(H I) + n(H II)
   κ_ν = α_ff(H⁻) * Pₑ * ------------------ * -----------------
                          n(H I) + n(H II)            ρ
This can be rewritten as: κ_ν = α_ff(H⁻) * Pₑ * n(H I) / ρ. We can rearrange this to get:
α_ff(H⁻) = κ_ν * ρ / (n(H I) * Pₑ).

I'm confident sure that the "more correct" version of the opacity equation should actually read 
κ_ν = α_ff(H⁻) * Pₑ * n(H I, n = 1) / ρ, and that the form provided by Gray (2005) implicitly
assumes that n(H I, n = 1) ~ n(H I). This seems to be a pretty good assumption since the partition 
function of H I is approximately equal to gₙ₌₁*exp(-Eₙ₌₁/(k*T)) = 2*exp(0) = 2, for most atmospheric
temperatures.
"""
function calc_H⁻_ff_absorption_coef(ν, T)
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

# I tabulated this data from Bell & Berrington (1987)
# each row corresponds to a separate λ (in Å)
_ff_λ_vals = [151890, 113918, 91134, 45567, 30378, 22784, 18227, 15189, 13019, 11392, 10126, 9113,
              7595, 6510, 5696, 5063, 4557, 3645, 3038, 2604, 2278, 1823]
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


if false # validate the accuracy of this table against figure 8.4 of Gray (2005)
    using Plots 
    labels = ["2520 K" "5040 K" "10080 K"]
    coef_vals = [view(_ff_table, :, 9), view(_ff_table, :, 4), view(_ff_table, :, 1)]
    p = plot(_ff_λ_vals, coef_vals, label = labels, xaxis=:log, yaxis=:log,
             legend=:topleft, xlim = [1e3,2e5], ylim = [1e-2, 500.],
             xlabel = "λ [Angstrom]", ylabel = "α_ff(H⁻) [1e-26 cm² per H atom per unit Pₑ]")
    display(p)
end


function check_Hminus_ff_values(target_precision = 0.01, verbose = true)
    # only look at the subtable that might be relevant
    sub_table = view(_ff_table, 2:20, 1:9)

    comp_λ = view(_ff_λ_vals, 2:20)
    comp_T = view(_ff_T_vals, 1:9)

    success = true
    ν_vals = (SSSynth.clight_cgs*1e8)./comp_λ
    for index = 1:length(comp_T)
        cur_T = comp_T[index]
        ref_coef = view(sub_table, :, index)
        calc_coef = calc_H⁻_ff_absorption_coef.(ν_vals, cur_T)
        precision = (abs.(calc_coef .- ref_coef)./ref_coef)
        max_err, max_err_ind = findmax(precision)
        if max_err > target_precision
            if verbose
                if success
                    println("There is a problem. The max error of polynomial should be: ",
                            target_precision*100, "%.")
                end
                println("Max error for T = ", cur_T, "K is: ", 100 * max_err,
                        " at λ = ", comp_λ[max_err_ind], " Å.")
            end
            success = false
        end
    end
    return success
end

@testset "H⁻ free-free opacity" begin
    @test check_Hminus_ff_values(0.0225)
end
