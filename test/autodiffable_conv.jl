using FiniteDiff, ForwardDiff

@testset "autodiffable_conv" begin
    # silly function involving convolution
    function f1(x)
        a = ones(5)*x[1]
        b = collect(1:7) * x[2]
        Korg.autodiffable_conv(a,b)
    end

    # silly function involving convolution where only one array contains duals
    function f2(x)
        a = ones(5)*x[1]*x[2]
        b = collect(1:7) 
        Korg.autodiffable_conv(a,b)
    end

    @testset for f in [f1, f2] 
        x = [2.0, 3.0]
        J_forward = ForwardDiff.jacobian(f,x) 
        J_finite = FiniteDiff.finite_difference_jacobian(f, x)
        @test J_forward ≈ J_finite rtol=1e-8
    end
end