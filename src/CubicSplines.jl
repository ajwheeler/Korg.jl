module CubicSplines
using LinearAlgebra

#=
This module contains modified functions from DataInterpolations.jl (license below).
We aren't depending on DataInterpolations because its dependencies are very heavy.
In the future, hopefully we can eliminate interpolation code and depend on a widely-used
well-tested and lightweight library.

Note, I've swapped the order of t and u (t is the abscissae/x-values, u are the y-values). The other 
major change is that I've simplified the types. We might at some point want to add analytic 
derivatives (DataInterpolations has these), but since ForwardDiff doesn't use ChainRules, there is 
no way to get them used in autodiff (until/unless you use a different library.)

Copyright (c) 2018: University of Maryland, Center for Translational Medicine.
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
=#

struct CubicSpline{tType,uType,hType,zType,T}
    t::tType
    u::uType
    h::hType
    z::zType
    extrapolate::Bool
    function CubicSpline(t, u, h, z, extrapolate)
        new{typeof(t),typeof(u),typeof(h),typeof(z),eltype(u)}(t, u, h, z, extrapolate)
    end
end

"""
    CubicSpline(xs, ys; extrapolate=false)

Construct a interpolant using `xs` and `ys` as the knot coordinates. Assumes `xs` is sorted. Apply
this object as a function to interpolate at any x value in the domain.
If `extrapolate` is false, x values outside [`xs[1]`, `xs[end]`] throw errors, if `extrapolate` is
true, the interpolant uses flat extrapolation, i.e. it returns the extreme value.
"""
function CubicSpline(t, u; extrapolate=false)
    n = length(t) - 1
    h = vcat(0, map(k -> t[k+1] - t[k], 1:length(t)-1), 0)
    dl = h[2:n+1]
    d_tmp = 2 .* (h[1:n+1] .+ h[2:n+2])
    du = h[2:n+1]
    tA = LinearAlgebra.Tridiagonal(dl, d_tmp, du)
    d = Vector{promote_type(eltype(u), eltype(h))}(undef, n + 1)
    d .= (i -> i == 1 || i == n + 1 ? 0 : 6(u[i+1] - u[i]) / h[i+1] - 6(u[i] - u[i-1]) / h[i]).(1:n+1)
    z = tA \ d
    CubicSpline(t, u, h[1:n+1], z, extrapolate)
end

function (A::CubicSpline{<:AbstractVector{<:Number}})(t::Number)
    if !(A.t[1] <= t <= A.t[end])
        if A.extrapolate #flat extrapolation
            if t < A.t[1]
                return A.u[1]
            else
                return A.u[end]
            end
        else
            throw(ArgumentError("Out-of-bounds value $(t) passed to interpolant. Must be between" *
                                " $(A.t[1]) and $(A.t[end])"))
        end
    end
    i = max(1, min(searchsortedlast(A.t, t), length(A.t) - 1))
    I = A.z[i] * (A.t[i+1] - t)^3 / (6A.h[i+1]) + A.z[i+1] * (t - A.t[i])^3 / (6A.h[i+1])
    C = (A.u[i+1] / A.h[i+1] - A.z[i+1] * A.h[i+1] / 6) * (t - A.t[i])
    D = (A.u[i] / A.h[i+1] - A.z[i] * A.h[i+1] / 6) * (A.t[i+1] - t)
    I + C + D
end

"""
    cumulative_integral!(out, A, t1, t2)

Given a curve described by the spine, A, Calculates the integral from t1 to t for all t = t1, t2, and
all spline knots in between.  So if `t1` is `A.t[1]` and `t2` is `A.t[end]`, `out` should have the
same length as `A.t`.
"""
function cumulative_integral!(out, A::CubicSpline, t1, t2)
    # the index less than or equal to t1
    idx1 = max(1, min(searchsortedlast(A.t, t1), length(A.t) - 1))
    # the index less than t2
    idx2 = max(2, min(searchsortedlast(A.t, t2), length(A.t) - 1))
    if A.t[idx2] == t2
        idx2 -= 1
    end

    out[1] = zero(eltype(A.u))
    for idx in idx1:idx2
        lt1 = idx == idx1 ? t1 : A.t[idx]
        lt2 = idx == idx2 ? t2 : A.t[idx+1]
        out[idx+1] = out[idx] + _integral(A, idx, lt2) - _integral(A, idx, lt1)
    end
end

function _integral(A::CubicSpline{<:AbstractVector{<:Number}}, idx::Number, t::Number)
    t1 = A.t[idx]
    t2 = A.t[idx+1]
    u1 = A.u[idx]
    u2 = A.u[idx+1]
    z1 = A.z[idx]
    z2 = A.z[idx+1]
    h2 = A.h[idx+1]
    (t^4 * (-z1 + z2) / (24 * h2) + t^3 * (-t1 * z2 + t2 * z1) / (6 * h2) +
     t^2 * (h2^2 * z1 - h2^2 * z2 + 3 * t1^2 * z2 - 3 * t2^2 * z1 - 6 * u1 + 6 * u2) / (12 * h2) +
     t * (h2^2 * t1 * z2 - h2^2 * t2 * z1 - t1^3 * z2 - 6 * t1 * u2 + t2^3 * z1 + 6 * t2 * u1) /
     (6 * h2))
end

end #module
