"""
This module provides a wrapper for `nlsolve`, from `NLsolve.jl`, which is a package for solving 
nonlinear systems of equations.

The wrapper extends `nlsolve` to allow the `residuals!` function to take additional fixed parameters
as a tuple of vectors.  This convieniently handles 
[julia's closure inneficiency problems](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured),
but more importantly facilitates the fast differentiation (with `ForwardDiff`) of the minimizer with 
respect to the parameters. 
"""
module NLSolver
using NLsolve, ForwardDiff

"""
The function wraps `NLsolve.nlsolve` to allow the `residuals!` function to take additional fixed 
parameters as a tuple of vectors. It has a specialized method to handle taking derivatives of the 
zero with respect to the parameters.

    
# arguments
- `residuals!`: Function that computes the residuals of the system of equations.
    The first argument is the vector of residuals, which is modified in place. The second argument 
    is the vector of variables, and the third argument is the vector of fixed parameters, which are
    fixed.
- `x0`: Initial guess for the solution.
- `p`: Fixed parameters.

All kwargs are passed to `NLsolve.nlsolve`.
    
Returns a pair whose first element is the zero (with dual number if appropriate) and whose 
second element is the `NLsolve.nlsolve` result (whose zero may not have dual numbers).
"""
function solve(residuals!, x0, p; nlsolve_kwargs...)
    r!(r, x) = residuals!(r, x, p)
    result = nlsolve(r!, x0; nlsolve_kwargs...)
    result.zero, result
end

function solve(residuals!, x0, ps::Vector{ForwardDiff.Dual{T, V, P}}; nlsolve_kwargs...) where {T, V, P}
    @info "specializing solve for ForwardDiff.Dual"
    @info "x0 type: $(typeof(x0))"

    p_values = [ForwardDiff.value.(p) for p in ps]
    p_partials = [ForwardDiff.partials.(p) for p in ps] 
    r!(r, x) = residuals!(r, x, p_values)
    sol = nlsolve(r!, x0; nlsolve_kwargs...)

    tmp = similar(sol.zero)
    # ∂r/∂x|x=sol.zero
    drdx = ForwardDiff.jacobian(tmp, sol.zero) do tmp, x
        residuals!(tmp, x, p_values)
    end
    # ∂r/∂p|x=sol.zero
    drdp = ForwardDiff.jacobian(tmp, p_values) do tmp, p
        residuals!(tmp, sol.zero, p)
    end
    dxdp = -(drdx \ drdp)
    partial_zero = dxdp * p_partials

    dual_zero = map(sol.zero, eachrow(partial_zero)) do v, p
        ForwardDiff.Dual{T}(v, p...)
    end

    dual_zero, sol
end
    
end