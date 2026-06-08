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

    function solve_system(p1, p2)
        if any(x -> x isa ForwardDiff.Dual, (p1, p2))
            vp1 = ForwardDiff.value(p1)
            vp2 = ForwardDiff.value(p2)
            res_val! = test_residuals_constructor(vp1, vp2)
            x_sol, conv, inf_norm = Korg.clipped_newton(res_val!, [2.0, 3.0])
            return Korg.implicit_jacobian_diff(test_residuals_constructor, x_sol, p1, p2)
        end
        res! = test_residuals_constructor(p1, p2)
        # start at [1.0, 1.0] to get the zero we want. (This system has two solutions.)
        x_sol, conv, inf_norm = Korg.clipped_newton(res!, [1.0, 1.0])
        return x_sol
    end

    p1_val = 8.0
    p2_val = 7.0

    # basic clipped_newton behavior (no duals)
    x_val = solve_system(p1_val, p2_val)
    @test x_val ≈ [2.0, 3.0]

    # differentiation wrt [p1 ; p2]
    J = ForwardDiff.jacobian(p -> solve_system(p[1], p[2]), [p1_val, p2_val])
    expected_J = [-1/7 2/7;
                  4/7 -1/7]

    @test J ≈ expected_J
end
