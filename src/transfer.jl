"""
    _ray_approximate_transfer_integral(τ, m, b)

The exact solution to \$\\int (m\\tau + b) \\exp(-\\tau)\$ d\\tau\$.
"""
function _ray_approximate_transfer_integral(τ, m, b)
    -exp(-τ) * (m*τ + b + m)
end

"""
    _plane_parallel_approximate_transfer_integral(τ, m, b)

The exact solution to \$\\int (m\\tau + b) E_2(\\tau)\$ d\\tau\$.

The exponential integral function, expint, captures the integral over the disk of the star to 
get the emergent astrophysical flux. You can verify it by substituting the variable of integration 
in the exponential integal, t, with mu=1/t.
"""
function _plane_parallel_approximate_transfer_integral(τ, m, b)
    1/6 * (τ*exponential_integral_2(τ)*(3b+2m*τ) - exp(-τ)*(3b + 2m*(τ+1)))
end

"""
    transfer_integral(τ, S; plane_parallel=true)

Compute exactly the solution to the transfer integral obtained be linearly interpolating the source 
function, `S` across optical depths `τ`, without approximating factor of exp(-τ) or E₂(τ).
"""
function transfer_integral(τ, S; plane_parallel=true)
    @assert size(τ) == size(S)
    indef_integral = if plane_parallel
        _plane_parallel_approximate_transfer_integral
    else
        _ray_approximate_transfer_integral
    end
    I = 0
    for i in 1:length(τ)-1
        m = (S[i+1] - S[i])/(τ[i+1] - τ[i])
        b = S[i] - m*τ[i]
        I += (indef_integral(τ[i+1], m, b) - indef_integral(τ[i], m, b))
    end
    I
end

"""
TODO
"""
function spherical_transfer(R, radii, α, S, μ_surface_grid)
    r0 = radii[end] #radii of inner and outer atmosphere layers
    
    I_type = typeof(promote(radii[1], α[1], S[1], μ_surface_grid[1])[1])
    I = Matrix{I_type}(undef, size(α, 2), length(μ_surface_grid))
    
    for (μ_ind, μ_surface) in enumerate(μ_surface_grid)
        # impact parameter of ray
        b = R * sqrt(1 - μ_surface^2)
        
        #calculate the index of the lowest layer the ray passes through
        i = argmin(abs.(radii .- b))
        if radii[i] < b
            i -= 1
        end
        if i == 1 || i == 0
            I[:, μ_ind] = zeros(size(α, 2))
            continue
        end #TODO eliminate for i == 1?
        
        #ds is the path length through each layer
        ds = -diff(@. sqrt([R ; radii[1:i]]^2 - b^2)) 
        dτ = α[1:i, :] .* ds
        τ = cumsum(dτ, dims=1)
        
        I[:, μ_ind] = [transfer_integral(τ[:, j], S[1:i, j]; plane_parallel=false)
                       for j in 1:size(τ, 2)]
        I[:, μ_ind] += if b > r0 
            #if the ray never leaves the model atmosphere, include 
            #the contribution from the other side of the star
            τ_prime = τ[end:end, :] .+ [cumsum(reverse(dτ[2:end, :], dims=1), dims=1) 
                                        ; zeros(size(τ[end:end, :]))]
            [transfer_integral(τ_prime[:, j], S[1:i, j]; plane_parallel=false)
             for j in 1:size(τ, 2)]
        else #otherwise assume I=S at atmosphere lower boundary
            @. exp(-τ[end, :]) * S[end, :]
        end
    end
    I
end

"""
    trapezoid_rule(xs, fs)

Approximate the integral from x₁ to x₂ of f(x) with the trapezoid rule given x-values `xs` and f(x)
values `fs`.

This should be good enough to numerically solve the transport equation, since model atmospheres
usually have carefully chosen knots.  We probably want to add higher-order aproximations later.
"""
function trapezoid_rule(xs, fs)
    Δs = diff(xs)
    weights = [0 ; Δs] + [Δs ; 0]
    sum(0.5 * weights .* fs)
end

"""
    exponential_integral_2(x)

Approximate second order exponential integral, E_2(x).  This stiches together several series 
expansions to get an approximation which is accurate within 1% for all `x`.
"""
function exponential_integral_2(x) 
    if x < 1.1
        _expint_small(x)
    elseif x < 2.5
        _expint_2(x)
    elseif x < 3.5
        _expint_3(x)
    elseif x < 4.5
        _expint_4(x)
    elseif x < 5.5
        _expint_5(x)
    elseif x < 6.5
        _expint_6(x)
    elseif x < 7.5
        _expint_7(x)
    elseif x < 9
        _expint_8(x)
    else
        _expint_large(x)
    end
end

function _expint_small(x) 
    #euler mascheroni constant
    ℇ = 0.57721566490153286060651209008240243104215933593992
    1 + ((log(x) + ℇ - 1) + (-0.5 + (0.08333333333333333 + (-0.013888888888888888 + 
                                                            0.0020833333333333333*x)*x)*x)*x)*x
end
function _expint_large(x)
    invx = 1/x
    exp(-x) * (1 + (-2 + (6 + (-24 + 120*invx)*invx)*invx)*invx)*invx
end
function _expint_2(x)
    x -= 2
    0.037534261820486914 + (-0.04890051070806112 + (0.033833820809153176 + (-0.016916910404576574 + 
                                          (0.007048712668573576 -0.0026785108140579598*x)*x)*x)*x)*x
end
function _expint_3(x)
    x -= 3
    0.010641925085272673   + (-0.013048381094197039   + (0.008297844727977323   + 
            (-0.003687930990212144   + (0.0013061422257001345  - 0.0003995258572729822*x)*x)*x)*x)*x
end
function _expint_4(x)
    x -= 4
    0.0031982292493385146  + (-0.0037793524098489054  + (0.0022894548610917728  + 
            (-0.0009539395254549051  + (0.00031003034577284415 - 8.466213288412284e-5*x )*x)*x)*x)*x
end
function _expint_5(x)
    x -= 5
    0.000996469042708825   + (-0.0011482955912753257  + (0.0006737946999085467  +
            (-0.00026951787996341863 + (8.310134632205409e-5   - 2.1202073223788938e-5*x)*x)*x)*x)*x
end
function _expint_6(x)
    x -= 6
    0.0003182574636904001  + (-0.0003600824521626587  + (0.00020656268138886323 + 
            (-8.032993165122457e-5   + (2.390771775334065e-5   - 5.8334831318151185e-6*x)*x)*x)*x)*x
end
function _expint_7(x)
    x -= 7
    0.00010350984428214624 + (-0.00011548173161033826 + (6.513442611103688e-5   + 
            (-2.4813114708966427e-5  + (7.200234178941151e-6   - 1.7027366981408086e-6*x)*x)*x)*x)*x
end
function _expint_8(x)
    x -= 8
    3.413764515111217e-5   + (-3.76656228439249e-5    + (2.096641424390699e-5   + 
            (-7.862405341465122e-6   + (2.2386015208338193e-6  - 5.173353514609864e-7*x )*x)*x)*x)*x
end

