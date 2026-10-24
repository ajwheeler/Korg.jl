"""
    SharedKnotsLinearInterpolants(knots, values)

A collection of linear interpolants that share the same knots. `values[j]` holds the
values of the `j`th function at `knots`. The knots must be strictly increasing, with at least two
points.

Evaluating every function at the same `x` with [`interpolate_all!`](@ref) finds the interval and
weight once, rather than once per function.
"""
struct SharedKnotsLinearInterpolants{K,V}
    knots::K
    values::V
end

"""
    interpolate_all!(out, itps::SharedKnotsLinearInterpolants, x)

Set `out[j]` to the `j`th interpolant in `itps` evaluated at `x`. Throws a `BoundsError` if `x` is
outside the knots (or `NaN`). Results are identical to those from
`Interpolations.linear_interpolation(knots, values[j])`, including at the knots.
"""
function interpolate_all!(out, itps::SharedKnotsLinearInterpolants, x)
    knots = itps.knots
    first(knots) <= x <= last(knots) || throw(BoundsError(knots, x))

    # Match the gridded interpolation's choice of interval at an exact knot.
    i = clamp(searchsortedfirst(knots, x, Base.Order.ForwardOrdering()) - 1,
              1, length(knots) - 1)
    weight = (x - knots[i]) / (knots[i+1] - knots[i])

    for j in eachindex(out, itps.values)
        vals = itps.values[j]
        out[j] = weight * vals[i+1] + (1 - weight) * vals[i]
    end
    return out
end
