import Dierckx

@testset "cubic spline interpolation" begin
    xs = 0:0.1:5
    ys = sin.(xs)

    ditp = Dierckx.Spline1D(xs, ys; k=3, bc="error", s=0.0)
    kitp = Korg.CubicSplines.CubicSpline(xs, ys)

    #these are only approximately the same because Dierckx is a BSpline library.
    x = 0:0.01:5
    @test Dierckx.evaluate(ditp, x)≈kitp.(x) rtol=1e-3

    @test_throws ArgumentError kitp(-1)
    @test_throws ArgumentError kitp(6)

    #test flat extrapolation
    itp = Korg.CubicSplines.CubicSpline(xs, ys; extrapolate=true)
    @test itp(-1) == 0
    @test itp(10) == sin(5)
end
