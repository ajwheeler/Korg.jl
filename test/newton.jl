@testset "Newton solver" begin
    function test_residuals_constructor(p1, p2)
        function residuals!(F, x)
            # p₁ = x₁ + 2x₂
            # p₂ = x₁² + x₂
            F[1] = x[1] + 2 * x[2] - p1
            F[2] = x[1]^2 + x[2] - p2
        end
        return residuals!
    end

    res! = test_residuals_constructor(8.0, 7.0)
    # start at [1.0, 1.0] to get the zero we want. (This system has two solutions.)
    x_sol, conv, inf_norm = Korg.clipped_newton(res!, [1.0, 1.0])
    @test conv
    @test x_sol ≈ [2.0, 3.0]
end
