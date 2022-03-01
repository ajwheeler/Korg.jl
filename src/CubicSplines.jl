module CubicSplines
using LinearAlgebra

#=
This module is contained modified functions from DataInterpolations.jl (license below).
We aren't depending on DataInterpolations because it's dependencies are very heavy.
In the future, hopefully we can eliminate interpolation code and depend on a widely-used
well-tested and lightweight library.

Note, I've swapped the order of t and u (t is the abscissae/x-values, u are the y-values). The other 
major change is that I've simplified the types. We might at some point want to add analytic 
derivatives (DataInterpolations has these), but since ForwardDiff doesn't use ChainRules, there is 
no way to get them used in autodiff (untill/unless you use a different library.)

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
    CubicSpline(t,u,h,z) = new{typeof(t),typeof(u),typeof(h),typeof(z),eltype(u)}(t,u,h,z)
end

"""
    CubicSpline(xs, ys)

Construct a interpolant using `xs` and `ys` as the knot coordinates. Assumes `xs` is sorted. Apply 
this object as a function to interpolate at any x value in the domain. If it is not within 
[`xs[1]`, `xs[end]`] throws an error.
"""
function CubicSpline(t,u)
    n = length(t) - 1
    h = vcat(0, map(k -> t[k+1] - t[k], 1:length(t)-1), 0)
    dl = h[2:n+1]
    d_tmp = 2 .* (h[1:n+1] .+ h[2:n+2])
    du = h[2:n+1]
    tA = LinearAlgebra.Tridiagonal(dl,d_tmp,du)
    d = map(i -> i == 1 || i == n + 1 ? 0 : 6(u[i+1] - u[i]) / h[i+1] - 6(u[i] - u[i-1]) / h[i], 1:n+1)
    z = tA\d
    CubicSpline(t,u,h[1:n+1],z)
end

function (A::CubicSpline{<:AbstractVector{<:Number}})(t::Number)
    if !(A.t[1] < t < A.t[end])
        throw(ArgumentError("Out-of-bounds value $(t) passed to interpolant. Must be between" * 
                            " $(A.t[1]) and $(A.t[end])")) 
    end
    i = max(1, min(searchsortedlast(A.t, t), length(A.t) - 1))
    I = A.z[i] * (A.t[i+1] - t)^3 / (6A.h[i+1]) + A.z[i+1] * (t - A.t[i])^3 / (6A.h[i+1])
    C = (A.u[i+1]/A.h[i+1] - A.z[i+1]*A.h[i+1]/6)*(t - A.t[i])
    D = (A.u[i]/A.h[i+1] - A.z[i]*A.h[i+1]/6)*(A.t[i+1] - t)
    I + C + D
end

end #module
