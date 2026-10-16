# tests of things in src/shared_knots_linear_interpolants.jl
using Interpolations: linear_interpolation

@testset "SharedKnotsLinearInterpolants" begin
    @testset "matches linear_interpolation" begin
        # Compare directly with one interpolator per function. Different values for each function
        # check that only the interval and weights are shared, not the values themselves.
        for T in (Float32, Float64)
            grids = (T[1, 1.25, 2, 4],
                     Korg.Wavelengths([(5000.0, 5001.0, 0.25), (5004.0, 5005.0, 0.5)]))
            for knots in grids
                vals = [T[1 + i / 10 for i in eachindex(knots)],
                    T[2 + 1 / i for i in eachindex(knots)]]
                itps = [linear_interpolation(knots, row) for row in vals]
                shared = Korg.SharedKnotsLinearInterpolants(knots, vals)
                # Nonuniform knots and separated windows exercise interval selection.
                # Include endpoints, midpoints, and adjacent floats to catch boundary errors.
                positions = [collect(knots); (knots[1:end-1] .+ knots[2:end]) ./ 2]
                append!(positions,
                        [x for k in knots
                         for x in (prevfloat(k), nextfloat(k))
                         if first(knots) <= x <= last(knots)])

                for x in positions
                    expected = [itp(x) for itp in itps]
                    out = similar(expected)
                    Korg.interpolate_all!(out, shared, x)
                    @test isequal(out, expected)
                end
                # Clamping the interval index must not silently allow extrapolation.
                for x in (prevfloat(first(knots)), nextfloat(last(knots)), NaN)
                    @test_throws BoundsError itps[1](x)
                    @test_throws BoundsError Korg.interpolate_all!(zeros(2), shared, x)
                end
            end
        end
    end

    @testset "automatic differentiation" begin
        # Compare values and derivatives with linear_interpolation. The two parameters are the
        # evaluation point and a scale applied to the values, respectively.
        knots = [1.0, 1.25, 2.0, 4.0]
        function evaluate(p, shared)
            vals = [p[2] .* [1.0, 3.0, 2.0, 6.0], p[2] .* [2.0, 4.0, 7.0, 5.0]]
            if shared
                Korg.interpolate_all!(similar(p, 2),
                                      Korg.SharedKnotsLinearInterpolants(knots, vals), p[1])
            else
                [linear_interpolation(knots, row)(p[1]) for row in vals]
            end
        end
        # At an interior knot, the selected interval determines the derivative wrt x.
        for x in (1.0, 1.1, 1.25, prevfloat(2.0), 2.0, nextfloat(2.0), 3.5)
            p = [x, 1.3]
            @test evaluate(p, true) == evaluate(p, false)
            @test ForwardDiff.jacobian(p -> evaluate(p, true), p) ≈
                  ForwardDiff.jacobian(p -> evaluate(p, false), p)
        end
    end
end
