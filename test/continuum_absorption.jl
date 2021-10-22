using Interpolations: bounds

include("absorption_comparison_funcs.jl")
include("utilities.jl")

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
    nH_I_div_partition = nH_I/2.0

    linear_abs_coef = Korg.ContinuumAbsorption.Hminus_ff([ν], T, nH_I_div_partition, ne)[1]
    Pₑ = ne * Korg.kboltz_cgs * T
    linear_abs_coef/(nH_I * Pₑ * 1e-26)
end

# the first axis of ref_val_matrix should have the same length as λ_vals_cm
# the second axis of ref_val_matrix should have the same length as T_vals
function _compare_against_table(ref_val_matrix, λ_vals_cm, T_vals, calc_func, target_precision,
                                verbose = true)
    @assert ((length(size(ref_val_matrix)) == 2) &
             (length(size(λ_vals_cm)) == 1) &
             (length(size(T_vals)) == 1))
    @assert size(ref_val_matrix)[1] == size(λ_vals_cm)[1]
    @assert size(ref_val_matrix)[2] == size(T_vals)[1]

    ν_vals = Korg.c_cgs ./ λ_vals_cm

    success = true
    for index = 1:length(T_vals)
        cur_T = T_vals[index]
        ref_vals = view(ref_val_matrix, :, index)
        calc_vals = calc_func.(ν_vals, cur_T)
        precision = (abs.(calc_vals .- ref_vals)./ref_vals)
        max_err, max_err_ind = findmax(precision)
        if max_err > target_precision
            if verbose
                if success
                    @info string("There is a problem. The max error of should be: ",
                                 target_precision*100, "%.")
                end
                @info string("Max error for T = ", cur_T, "K is: ", 100 * max_err,
                             "% at λ = ", (λ_vals_cm[max_err_ind] * 1.e8), " Å.")
            end
            success = false
        end
    end
    return success

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
    _compare_against_table(view(_ff_table, 2:20, 1:9), view(_ff_λ_vals, 2:20) ./ 1e8,
                           view(_ff_T_vals, 1:9), _bell_berrington_87_Hminus_ff_absorption_coef,
                           target_precision, verbose)
end

@testset "H⁻ free-free absorption" begin
    @test check_Hminus_ff_values(0.0225)
    # we only have measurements for b and c
    @testset "Gray (2005) Fig 8.5$panel comparison" for panel in ["b", "c"]
        calculated, ref = Gray_opac_compare.Gray05_comparison_vals(panel,"Hminus_ff")
        @test all(calculated .≥ 0.0)
        @test all(abs.(calculated - ref) .≤ Gray_opac_compare.Gray05_atols[panel])
    end
end

"""
    _Hminus_bf_cross_section_Gray(λ)

Compute the H⁻ bound-free cross-section at a given wavelength (specified in Å). The cross-section 
has units of megabarns per H⁻ particle and does NOT include a correction for stimulated emission.

This function uses the polynomial provided in equation 8.11 of Gray (2005), that fits the tabulated
data from Wishart (1979). While Gray (2005) claims that the polynomial fits the data with 0.2%
precision for 2250 Å ≤ λ ≤ 15000 Å, in practice we find that it fits the data to better than 0.25%
precision.
"""
function _Hminus_bf_cross_section_Gray(λ::Real)
    # we need to somehow factor out this bounds checking
    if !(2250 <= λ <= 15000.0)
        throw(DomainError(λ, "The wavelength must lie in the interval [2250 Å, 15000 Å]"))
    end

    λ2 = λ*λ
    λ3 = λ*λ2
    λ4 = λ*λ3
    λ5 = λ*λ4
    λ6 = λ*λ5

    αbf_H⁻ = (1.99654 - 1.18267e-5 * λ + 2.64243e-6 * λ2 - 4.40524e-10 * λ3 + 3.23992e-14 * λ4
              - 1.39568e-18 * λ5 + 2.78701e-23 * λ6)
    αbf_H⁻
end

function check_Hminus_bf_values(target_precision = 0.002, verbose = true)
    # compare the tabulated data from Wishart (1979) against the polynomial fits from Gray (2005)
    
    _bf_table = Korg.ContinuumAbsorption._Hminus_bf_table
    _bf_λ_vals = view(_bf_table, :, 1)
    _bf_α_vals = view(_bf_table, :, 2)

    comp_λ_vals = view(_bf_λ_vals, 5:56)
    comp_α_vals = view(_bf_α_vals, 5:56)
    @assert comp_λ_vals[1] == 2250.0
    @assert comp_λ_vals[end] == 15000.0

    coefs = _Hminus_bf_cross_section_Gray.(comp_λ_vals)
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


@testset "H⁻ bound-free absorption" begin
    @test check_Hminus_bf_values(0.0025)
    @testset "Gray (2005) Fig 8.5$panel comparison" for panel in ["a", "b", "c"]
        calculated, ref = Gray_opac_compare.Gray05_comparison_vals(panel,"Hminus_bf")
        @test all(calculated .≥ 0.0)
        @test all(abs.(calculated - ref) .≤ Gray_opac_compare.Gray05_atols[panel])
    end
    @testset "Extreme wavelengths" begin
        # this tests the absorption function at wavelengths outside of the Wishart (1979) table

        # choose arbitrary physical values:
        nH_I = 3.0e16
        nH_I_div_partition = nH_I / 2.0
        ne = nH_I / 100.0
        T = 7800.0

        # determine the minimum and maximum λs in the table:
        min_tabulated_λ = Korg.ContinuumAbsorption._Hminus_bf_table[1,1]
        max_tabulated_λ = Korg.ContinuumAbsorption._Hminus_bf_table[end, 1]
        # determine the ionization λ for H⁻:
        ion_energy = Korg.ContinuumAbsorption._H⁻_ion_energy
        Å_per_eV = 1e8 * (Korg.hplanck_eV * Korg.c_cgs)
        max_λ_ionize = Å_per_eV/ion_energy

        # now we are ready for the tests:

        # first, check that a bounds error is thrown below min_tabulated_λ
        @test_throws DomainError Korg.ContinuumAbsorption.Hminus_bf(
            [Korg.c_cgs/(1e-8*0.5*min_tabulated_λ)], T, nH_I_div_partition, ne, ion_energy;
            error_oobounds = true
        )[1]

        # next, check that the linear absorption coefficient between max_tabulated_λ and
        # max_λ_ionize is between the absorption coefficient at max_tabulated_λ and 0
        α_max_tabulated_λ = Korg.ContinuumAbsorption.Hminus_bf(
            [Korg.c_cgs/(1e-8*max_tabulated_λ)], T, nH_I_div_partition, ne, ion_energy
        )[1]

        α_test = Korg.ContinuumAbsorption.Hminus_bf(
            [Korg.c_cgs/(1e-8*0.5*(max_tabulated_λ+max_λ_ionize))], T, nH_I_div_partition, ne,
            ion_energy
        )[1]
        @test (α_max_tabulated_λ > α_test) && (α_test > 0.0)

        # finally, check that the linear absorption coefficient at λ > max_λ_ionize is zero
        @test 0.0 == Korg.ContinuumAbsorption.Hminus_bf(
            [Korg.c_cgs/(1e-8*2*max_λ_ionize)], T, nH_I_div_partition, ne, ion_energy
        )[1]
    end
end


# gives the atomic absorption coefficient (or cross section)
function check_Heminus_ff_absorption(target_precision, verbose = true)

    # taken from Table 2 of John (1994):
    _Heminus_table =
        [27.979 30.907 35.822 39.921 43.488 46.678 49.583 52.262 54.757 63.395 70.580;
         15.739 17.386 20.151 22.456 24.461 26.252 27.882 29.384 30.782 35.606 39.598;
         10.074 11.128 12.897 14.372 15.653 16.798 17.838 18.795 19.685 22.747 25.268;
          4.479  4.947  5.733  6.387  6.955  7.460  7.918  8.338  8.728 10.059 11.147;
          2.520  2.784  3.226  3.593  3.910  4.193  4.448  4.681  4.897  5.632  6.234;
          1.614  1.783  2.065  2.299  2.502  2.681  2.842  2.990  3.126  3.592  3.979;
          1.121  1.239  1.435  1.597  1.737  1.860  1.971  2.073  2.167  2.490  2.765;
          0.632  0.698  0.808  0.899  0.977  1.045  1.108  1.165  1.218  1.405  1.574;
          0.405  0.447  0.518  0.576  0.625  0.670  0.710  0.747  0.782  0.910  1.030;
          0.282  0.311  0.360  0.400  0.435  0.466  0.495  0.522  0.547  0.643  0.737;
          0.159  0.176  0.204  0.227  0.247  0.266  0.283  0.300  0.316  0.380  0.444;
          0.102  0.113  0.131  0.147  0.160  0.173  0.186  0.198  0.210  0.258  0.305;
          0.072  0.079  0.092  0.103  0.114  0.124  0.133  0.143  0.152  0.190  0.227;
          0.053  0.059  0.069  0.077  0.086  0.094  0.102  0.109  0.117  0.148  0.178;
          0.041  0.045  0.053  0.061  0.067  0.074  0.081  0.087  0.094  0.120  0.145;
          0.033  0.036  0.043  0.049  0.055  0.061  0.066  0.072  0.078  0.100  0.121] .* 1e-26

    _θ_vals_He⁻_ff_john94 = [0.5 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.8 3.6]
    _λ_vals_He⁻_ff_john94 = [15.1878, 11.3909, 9.1127, 6.0751, 4.5564, 3.6451, 3.0376, 2.2782,
                             1.8225,  1.5188, 1.1391, 0.9113, 0.7594, 0.6509, 0.5695, 0.5063] #μm

    T_vals = 5040.0 ./ _θ_vals_He⁻_ff_john94
    ν_vals = (Korg.c_cgs*1e4) ./_λ_vals_He⁻_ff_john94

    # For the sake of easy comparisons, let's pick the following semi-realistic values:
    nₑ = 1e13 # cm⁻³
    nHe_I = 8.5e13 # cm⁻³
    # for simplicity, assume that the He I partition function is always just equal to 1 and (thus
    # the number density in the ground state is equal to the total number density of He I)
    nHe_I_div_U, nHe_I_gs = nHe_I, nHe_I

    # this includes the correction for stimulated emission
    Pₑ = nₑ*Korg.kboltz_cgs.*T_vals
    ref_linear_absorption_vals = _Heminus_table .* nHe_I .* Pₑ

    # the extra handling in ν is necessary due to roundoff errors.
    calc_func(ν,T) = Korg.ContinuumAbsorption.Heminus_ff(
        [ifelse(ν == 5.921241516887221e14, prevfloat(ν), ν)], T, nHe_I_div_U, nₑ
    )[1]
    _compare_against_table(view(ref_linear_absorption_vals, :, 1:9), _λ_vals_He⁻_ff_john94 ./ 1e4,
                           view(T_vals, 1:9), calc_func, target_precision, verbose)
end

@testset "He⁻ free-free absorption" begin
    @test check_Heminus_ff_absorption(0.031)
    # this really only amounts to a sanity check because the absolute tolerance is of the same
    # magnitude as the actual values
    @testset "Gray (2005) Fig 8.5$panel comparison" for panel in ["b", "c"]
        calculated, ref = Gray_opac_compare.Gray05_comparison_vals(panel,"Heminus_ff")
        @test all(calculated .≥ 0.0)
        @test all(abs.(calculated - ref) .≤ Gray_opac_compare.Gray05_atols[panel])
    end
end

function check_H2plus_ff_and_bf_absorption(target_precision, verbose = true)
    # Table II from Bates (1952) gives the "absorption coefficient" (corrected for stimulated
    # emission) in units of 1e-39 cm⁻¹ / (H atom/cm³) / (H⁺ ion/cm³).
    #
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
    _T_vals = [2.5e3, 3.0e3, 3.5e3, 4.0e3, 5.0e3, 6.0e3, 7.0e3, 8.0e3, 1.0e4, 1.2e4]
    # the rows are different wavenumbers (in cm⁻¹)
    _wavenumbers = [ 4000,  5000,  6000,  7000,  8000,  9000, 10000, 12000, 14000, 16000, 18000,
                    20000, 22000, 24000, 26000]

    # to match the implicit assumption in Bates (1952) that ndens(H I) = ndens(H I, n = 1)
    partition_val = 2.0

    # pick stupid dummy values (we're just going to remove them later):
    nH_I = 1e15
    nH_II = 1e13

    nH_I_divU = nH_I/partition_val
    const_factor = 1e39 /(nH_I * nH_II)
    calc_func(ν,T) =  const_factor * Korg.ContinuumAbsorption.H2plus_bf_and_ff([ν], T, nH_I_divU,
                                                                                nH_II)[1]

    λ_cgs = 1.0 ./_wavenumbers
    _compare_against_table(_table, λ_cgs, _T_vals, calc_func, target_precision, verbose)
end

@testset "combined H₂⁺ ff and bf absorption" begin
    @test check_H2plus_ff_and_bf_absorption(0.015)
    @testset "Gray (2005) Fig 8.5$panel comparison" for panel in ["b"]
        calculated, ref = Gray_opac_compare.Gray05_comparison_vals(panel,"H2plus")
        @test all(calculated .≥ 0.0)
        @test all(abs.(calculated - ref) .≤ Gray_opac_compare.Gray05_atols[panel])
    end
end

"""
    compare_gauntff_kurucz([rtol])

Compare the free-free gaunt factor to the reference values tabulated in section 5.1 of Kurucz
(1970). These reference tabulated values are less accurate than the values we're actually using
(since these values were derived from a figure in Karsas and Latter (1961)), but they're adequate
for testing purposes.
"""
function compare_gauntff_kurucz(rtol = 0.15)

    # First, get bounds of log₁₀(u) = log₁₀(RydbergH*Z²/(k*T)) & log₁₀(γ²) = log₁₀(Rydberg*Z²/(k*T)),
    # in which our actual data is known.
    _log10_u_bounds, _log10_γ2_bounds = bounds(Korg.ContinuumAbsorption._gauntff_interpolator.itp)
    log10_u_bounds = Korg.Interval(_log10_u_bounds...)
    log10_γ2_bounds = Korg.Interval(_log10_γ2_bounds...)


    # Next, initialize the reference data
    log10_γ2 = [-3.0, -2.5, -2.0, -1.5, -1.0, -0.5,  0.0,  0.5,  1.0,  1.5,  2.0]
    log10_u = [-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5,  0.0,  0.5,  1.0,  1.5]
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
    assert_allclose_grid(
        Korg.ContinuumAbsorption._gauntff_interpolator.(cmp_log10_u, cmp_log10_γ2'),
        ref_gaunt_ff[ref_u_slc, ref_γ2_slc],
        [("log₁₀u", cmp_log10_u), ("log₁₀γ²", cmp_log10_γ2)];
        rtol = rtol, atol = 0.0
    )
end

"""
    gray_H_I_ff_absorption_coef(λ, T, [gaunt_func])

Computes the H I free-free linear absorption coefficient per neutral Hydrogen atom by adapting 
equation 8.10 from Gray (2005). This includes the correction for stimulated emission under LTE.
In contrast to `Korg.ContinuumAbsorption.H_I_ff`, this implementation directly combines the Saha 
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
function gray_H_I_ff_absorption_coef(λ, T, gaunt_func = nothing)

    # adapt equation 8.10
    R = 1.0968e5 # cm⁻¹
    I = Korg.hplanck_eV * Korg.c_cgs * R # eV
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
        β = 1.0/(Korg.kboltz_eV * T)
        ν = Korg.c_cgs * 1e8 / λ
        log_u = log10(Korg.hplanck_eV * ν * β)
        log_γ2 = log10(Korg.RydbergH_eV * Z2 * β)
        gaunt_func(log_u, log_γ2)
    end

    uncorrected = α₀ *λ^3 * gaunt_ff * loge * 10.0^(-θ*I)/ (2*θ*I)
    # according to equation 8.18, we need to correct the value for stimulated emission
    uncorrected*(1 - 10.0^(-χ_λ * θ))
end


@testset "H I free-free absorption" begin
    @testset "Kurucz (1970) Free-Free Gaunt Factor comparison" begin
        @test compare_gauntff_kurucz(0.15)
    end

    @testset "Gray (2005) implementation comparison" begin
        # Compare the different formulations of the free-free absorption under the conditions of
        # the panels of figure 8.5 of Gray (2005). Outside of these conditions, it's unclear where
        # Gray's approximation for the gaunt factor breaks down.
        # This comparison is much more constraining than directly comparing against the relevant
        # curve extracted from the panel because the relevant curve includes both free-free and
        # bound-free absorption. It's also much more constraining when H I ff absorption is
        # subdominant
        χs = [(Korg.RydbergH_eV, -1.0, -1.0)]
        # implicitly assumed by the Gray implementation
        Us = Dict([Korg.literals.H_I=>(T -> 2.0), Korg.literals.H_II=>(T -> 1.0)])

        λ_vals = [3e3 + (i-1)*250 for i = 1:69] # equally spaced vals from 3e3 Å through 2e4 Å
        ν_vals = Korg.c_cgs*1e8./λ_vals

        for (logPₑ, T) in [(1.08, 5143.0), (1.77, 6429.0), (2.50, 7715.0), (2.76, 11572.0)]
            ref_absorption_coef = gray_H_I_ff_absorption_coef.(λ_vals, T)

            nₑ = 10.0^logPₑ / (Korg.kboltz_cgs * T)
            nH_total = nₑ * 100.0 # this is totally arbitrary

            wII, wIII = Korg.saha_ion_weights(T, nₑ, 1, χs, Us)
            nH_I = nH_total / (1 + wII) 
            nH_II = nH_total * wII / (1 + wII)
            # recall that H I ff actually refers to: photon + e⁻ + H II -> e⁻ + H II
            absorption_coef = Korg.ContinuumAbsorption.H_I_ff(ν_vals, T, nH_II, nₑ) / nH_I
            @test assert_allclose_grid(absorption_coef, ref_absorption_coef, [("λ", λ_vals, "Å")];
                                       rtol = 0.032, atol = 0)
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
This implicitly assumes that the H I partition function is always equal to 2.
"""
function gray_H_I_bf_absorption_coef(λ, T)
    # equation 8.8 gives the formula without stimulated emission corrections
    # λ must have units of Å and T must have units of K

    θ = 5040/T # eV^-1
    photon_energy_eV = Korg.hplanck_eV * Korg.c_cgs * 1e8/λ

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

@testset "H I bound-free absorption" begin
    @testset "Gray (2005) implementation comparison" begin
        # Compare the different formulations of the bound-free absorption under the conditions of
        # the panels of figure 8.5 of Gray (2005). Outside of these conditions, it's unclear
        # if/where Gray's approximation for the gaunt factor breaks down.
        # This comparison is much more constraining than directly comparing against the relevant
        # curve extracted from the panel because the relevant curve includes both free-free and
        # bound-free absorption. It's also much more constraining when H I bf absorption is
        # subdominant

        H_I_ion_energy = Korg.RydbergH_eV # use this to decouple test from ionization energies
        H_I_partition_val = 2.0 # implicitly assumed by the Gray implementation

        λ_vals = [3e3 + (i-1)*250 for i = 1:69] # equally spaced vals from 3e3 Å through 2e4 Å
        ν_vals = Korg.c_cgs*1e8./λ_vals

        for (logPₑ, T) in [(1.08, 5143.0), (1.77, 6429.0), (2.50, 7715.0), (2.76, 11572.0)]
            ref_absorption_coef = gray_H_I_bf_absorption_coef.(λ_vals, T)

            nₑ = 10.0^logPₑ / (Korg.kboltz_cgs * T)
            nH_I = nₑ * 100.0 # this is totally arbitrary

            absorption_coef = Korg.ContinuumAbsorption.H_I_bf(ν_vals, T, nH_I/H_I_partition_val,
                                                              H_I_ion_energy) / nH_I
            @test maximum(abs.(absorption_coef - ref_absorption_coef)/absorption_coef) < 0.004
            @test all(absorption_coef .≥ 0.0)
        end
    end

    @testset "Integral-Summation Equivalence" begin
        # Korg.ContinuumAbsorption.H_I_bf approximates the sum of the absorption contributions from
        # bound-free transitions originating from high energy levels with an integral. These tests
        # check that the integral does a reasonable job reproducing the direct sum.

        H_I_ion_energy = Korg.RydbergH_eV # use this to decouple test from ionization energies
        H_I_partition_val = 2.0 # implicitly assumed by the Gray implementation

        λ_vals = [3e3 + (i-1)*500 for i = 1:35] # equally spaced vals from 3e3 Å through 2e4 Å
        ν_vals = Korg.c_cgs*1e8./λ_vals

        nstop_sum = 100

        for (logPₑ, T) in [(1.08, 5143.0), (1.77, 6429.0), (2.50, 7715.0), (2.76, 11572.0)]
            ref_absorption_coef = gray_H_I_bf_absorption_coef.(λ_vals, T)

            nₑ = 10.0^logPₑ / (Korg.kboltz_cgs * T)
            nH_I = nₑ * 100.0 # this is totally arbitrary

            α_integral = Korg.ContinuumAbsorption.H_I_bf(ν_vals, T, nH_I/H_I_partition_val,
                                                         H_I_ion_energy)
            α_sum = Korg.ContinuumAbsorption.H_I_bf(ν_vals, T, nH_I/H_I_partition_val,
                                                    H_I_ion_energy, nstop_sum, false)
            #println(maximum(abs.(α_integral - α_sum)/α_sum))
            @test maximum(abs.(α_integral - α_sum)/α_sum) < 0.003
        end
    end
end


# this is only included here for completeness.
# Given the strong consistency when comparing the bf and ff values individually against the
# alternative implementations described in Gray (2005), I'm unconcerned by the need for large
# tolerances.
@testset "combined H I bf and ff absorption" begin
    @testset "Gray (2005) Fig 8.5$panel comparison" for (panel, atol) in [("b",0.035),
                                                                          ("c",0.125),
                                                                          ("d",38.5)]
        calculated, ref = Gray_opac_compare.Gray05_comparison_vals(panel,"H")
        @test all(calculated .≥ 0.0)
        @test assert_allclose(calculated, ref; rtol = 0.0, atol = atol)
    end
end


@testset "TOPbase bound-free absorption" begin
    T = 7800.0 #K, this is fairly arbitrary
    ndens_species = 3.0e16 #cm⁻³, this is fairly arbitrary

    @testset "$species_name comparison" for species_name in ["H_I", "He_II"]
        λ_vals = OP_compare._dflt_λ_vals(species_name)

        # compute the absorption coefficients using the TOPbase data
        hydrogenic_α_OP = OP_compare.calc_hydrogenic_bf_absorption_coef(λ_vals, T, ndens_species,
                                                                        species_name;
                                                                        use_OP_data = true)
        # compute the absorption coefficients using our function for hydrogenic atoms
        hydrogenic_α_dflt = OP_compare.calc_hydrogenic_bf_absorption_coef(λ_vals, T, ndens_species,
                                                                          species_name;
                                                                          use_OP_data = false)

        λ_comp_intervals = OP_compare._λ_comp_intervals(species_name)
        comp_ind = map(λ_vals) do λ
            any(λ_comp_intervals[:, 1] .<= λ .<= λ_comp_intervals[:, 2])
        end
        @test assert_allclose_grid(hydrogenic_α_OP[comp_ind], hydrogenic_α_dflt[comp_ind],
                                   [("λ", λ_vals[comp_ind], "Å"),];
                                   rtol = OP_compare._hydrogenic_rtol(species_name), atol = 0.0,
                                   err_msg = ("\n$(species_name) bf absorption coefficients " *
                                              "computed using data from the opacity project are " *
                                              "inconsistent with the results computed for a " *
                                              "hydrogenic atom"))

    end
end
