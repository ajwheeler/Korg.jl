@testset "moog style transfer" begin

using SpecialFunctions: expint
@testset "transfer" begin
    xs = 2:0.01:8
    @test Korg.MoogStyleTransfer.exponential_integral_2.(xs) ≈ expint.(2, xs) rtol=1e-3
end

@testset "trapezoid rule" begin
    #gaussian PDF should integral to 1.
    pdf(x) = exp(-1/2 * x^2) / sqrt(2π)
    xs = -10:0.1:10
    fs = pdf.(xs)

    partial_ints = similar(xs)
    Korg.MoogStyleTransfer.cumulative_trapezoid_rule!(partial_ints, xs, fs)

    naive_partial_integrals = map(1:length(xs)) do i
        Korg.MoogStyleTransfer.trapezoid_rule(xs[1:i], fs[1:i])
    end

    #do they match the analytic solution?
    @test Korg.MoogStyleTransfer.trapezoid_rule(xs, fs) ≈ 1.0  atol=1e-5
    @test partial_ints[end] ≈ 1.0 atol=1e-5

    #do they match each other?
    @test naive_partial_integrals ≈ partial_ints
end

@testset "synthesize spectrum" begin
    atm = read_model_atmosphere("data/sun.mod")
    legacy_sol = synthesize(atm, [], 5000, 5001; use_legacy_radiative_transfer=true, n_mu_points=50)
    sol = synthesize(atm, [], 5000, 5001 )
    @test assert_allclose_grid(legacy_sol.flux, sol.flux, [("λ" , sol.wavelengths, "Å")]; rtol=0.02)
end

end