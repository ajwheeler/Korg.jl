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