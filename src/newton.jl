# Newton solver with step clipping. Returns (x, converged, final_residuals_inf_norm).
function clipped_newton(residuals!, x0; tol=1e-8, max_iter=200, max_step=1.0)
    x = copy(x0)
    F = zeros(eltype(x), length(x))
    Fwork = similar(F)
    inf_norm = oftype(real(zero(eltype(x))), Inf) # TODO could this be streamlined?
    for _ in 1:max_iter
        residuals!(F, x)
        inf_norm = maximum(abs, F)
        if !isfinite(inf_norm)
            return x, false, inf_norm
        end
        if inf_norm < tol
            return x, true, inf_norm
        end

        # else, take a step

        # compute naïve Newton step
        # TODO would be nice to not allocate
        J = ForwardDiff.jacobian((Fout, xin) -> residuals!(Fout, xin), Fwork, x)
        all(isfinite, J) || return x, false, inf_norm
        step = try
            -(J \ F)
        catch
            return x, false, inf_norm
        end
        all(isfinite, step) || return x, false, inf_norm

        # clip the step so that the max diff in any dimension does not exceed max_step
        smax = maximum(abs, step)
        α = smax > max_step ? max_step / smax : one(eltype(step))
        x .+= α .* step
    end
    return x, false, inf_norm
end

"""
    implicit_jacobian_diff(residuals_constructor, x_sol, p...)

Differentiate the solution `x_sol` of an implicit system `R(x, p) = 0` with respect to the parameters
`p...` using the implicit function theorem.

  - `residuals_constructor(p...) -> residuals!(F, x)`: returns the in-place residuals function for the
    given parameter values.
  - `x_sol`: the converged Float64 solution.
  - `p...`: parameters of the system; at least one must be a `ForwardDiff.Dual`.
"""
function implicit_jacobian_diff(residuals_constructor, x_sol, p...)
    # get the primal parts of p...
    vp = map(x -> x isa AbstractArray ? ForwardDiff.value.(x) : ForwardDiff.value(x), p)

    # evaluate the residuals at x_sol with dual parameters to get the partials: ∂F = ∂R/∂p · ∂p
    F_dual = zeros(promote_type(eltype(x_sol), map(eltype, p)...), length(x_sol))
    residuals_constructor(p...)(F_dual, x_sol)
    ∂F = stack(i -> ForwardDiff.partials(F_dual[i]).values, eachindex(F_dual))'

    # Compute ∂R/∂x at the solution using Float64 values.
    tmp = similar(x_sol)
    drdx = ForwardDiff.jacobian((tmp, x) -> residuals_constructor(vp...)(tmp, x), tmp, x_sol)

    # Solve ∂x/∂p · ∂p = −(∂R/∂x)⁻¹ (∂R/∂p · ∂p) via the partials already embedded in F_dual.
    ∂x = -(drdx \ ∂F)

    T = ForwardDiff.tagtype(eltype(F_dual))
    [ForwardDiff.Dual{T}(x_sol[i], ∂x[i, :]...) for i in eachindex(x_sol)]
end