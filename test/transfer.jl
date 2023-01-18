@testset "transfer" begin

using SpecialFunctions: expint
using Korg.RadiativeTransfer: MoogStyleTransfer

@testset "transfer" begin
    xs = 2:0.01:8
    @test MoogStyleTransfer.exponential_integral_2.(xs) ≈ expint.(2, xs) rtol=1e-3
end

@testset "trapezoid rule" begin
    #gaussian PDF should integral to 1.
    pdf(x) = exp(-1/2 * x^2) / sqrt(2π)
    xs = -10:0.1:10
    fs = pdf.(xs)

    partial_ints = similar(xs)
    MoogStyleTransfer.cumulative_trapezoid_rule!(partial_ints, xs, fs)
    @test partial_ints[end] ≈ 1.0 atol=1e-5
end

@testset "synthesize spectrum with bezier transfer" begin
    atm = read_model_atmosphere("data/sun.mod")
    bezier_sol = synthesize(atm, [], format_A_X(), 5000, 5001; bezier_radiative_transfer=true, n_mu_points=50)
    sol = synthesize(atm, [], format_A_X(), 5000, 5001 )
    @test assert_allclose_grid(bezier_sol.flux, sol.flux, [("λ" , sol.wavelengths, "Å")]; rtol=0.02)
end

end