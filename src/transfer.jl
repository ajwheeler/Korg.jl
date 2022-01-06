"""
    radiative_transfer(atm::ModelAtmosphere, α, S, α_ref, mu_grid)

Returns the astrophysical flux at each wavelength.

inputs:
- `atm`: the model atmosphere.
- `α`: a matrix (atmospheric layers × wavelengths) containing the absoprtion coefficient
- `S`: the source fuction as a matrix of the same shape.
- `α_ref`: the continuum absoprtion coefficient at the reference wavelength (5000 Å). Used to 
   rescale the total absorption to match the model atmosphere. This value should be calculated by 
   Korg.
- `mu_grid`: (required if atm is a [`ShellAtmosphere`](@ref)) the values of μ at which to calculate 
   the surface intensity, which is integrated to obtain the astrophysical flux.
"""
function radiative_transfer(atm::PlanarAtmosphere, α, S, α_ref, mu_grid=nothing)
    τ5 = [l.tau_5000 for l in atm.layers] #τ at 5000 Å according to model atmosphere
    planar_transfer(α, S, τ5, α_ref)
end
function radiative_transfer(atm::ShellAtmosphere, α, S, α_ref, mu_grid)
    τ5 = [l.tau_5000 for l in atm.layers] #τ at 5000 Å according to model atmosphere
    radii = [atm.R + l.z for l in atm.layers]
    photosphere_correction = radii[1]^2 / atm.R^2
    #discard I, take F only
    photosphere_correction * spherical_transfer(α, S, τ5, α_ref, radii, mu_grid)[2]
end


"""
    planar_transfer(α, S, τ_ref, α_ref)

Returns the astrophysical flux. See [`radiative_transfer`](@ref) for an explantion of the arguments.
"""
function planar_transfer(α, S, τ_ref, α_ref)
    #this could possibly be sped up with explicit allocations and looping
    map(zip(eachcol(α), eachcol(S))) do (α_λ, S_λ)
        τ = cumulative_trapezoid_rule(log.(τ_ref), τ_ref .* α_λ ./ α_ref) 
        2π * transfer_integral(τ, S_λ; plane_parallel=true)
    end
end

"""
    spherical_transfer(α, S, τ_ref, α_ref, radii, μ_surface_grid)

Perform radiative transfer along rays emerging at the μ values in `μ_surface_grid` in a spherically
symmetric atmosphere resolved at radii `radii` [cm]. See [`radiative_transfer`](@ref) for an 
explantion of the arguments. Note that `radii` should be in decreasing order.

Returns `(flux, intensity)`, where `flux` is the astrophysical flux, and `intensity`, a matrix of 
shape (wavelengths × mu values), is the surface intensity as a function of μ.
"""
function spherical_transfer(α, S, τ_ref, α_ref, radii, μ_surface_grid)
    R, r0 = radii[1], radii[end] #lower bound of atmosphere

    #I is intensity as a funciton of μ and λ, filled in the loop below
    I_type = typeof(promote(radii[1], α[1], S[1], μ_surface_grid[1])[1])
    I = Matrix{I_type}(undef, size(α, 2), length(μ_surface_grid)) 

    #do radiative transfer along each ray
    for (μ_ind, μ_surface) in enumerate(μ_surface_grid)
        # impact parameter of ray
        b = R * sqrt(1 - μ_surface^2)
        
        #calculate the index of the lowest layer the ray passes through.
        #doing this with `findfirst` is messier at first and last index
        i = argmin(abs.(radii .- b)) 
        if radii[i] < b
            i -= 1
        end
        if i == 0 
            I[:, μ_ind] .= 0
            continue
        end 

        #geometric path-length correction
        ds_dr = @. radii[1:i] ./ sqrt(radii[1:i]^2 - b^2) 

        for j in 1:size(α, 2) #iterate over wavelength
            integrand = ds_dr .* τ_ref[1:i] .* α[1:i, j] ./ α_ref[1:i] 
            τ_λ = cumulative_trapezoid_rule(log.(τ_ref[1:i]), integrand)
            I[j, μ_ind] = transfer_integral(τ_λ, S[1:i, j]; plane_parallel=false)

            if b > r0
                #if the ray never leaves the model atmosphere, include the contribution from the 
                #other side of the star. This should probably by audited for off-by-one errors.
                τ_prime = τ_λ[end] .+ [cumsum(reverse(diff(τ_λ))) ; 0]
                I[j, μ_ind] += transfer_integral(τ_prime, S[1:i, j]; plane_parallel=false)
            else #otherwise assume I=S at atmosphere lower boundary.  This is a _tiny_ effect.
                I[j, μ_ind] += exp(-τ_λ[end]) * S[end, j]
            end
        end
    end
    #calculate 2π∫μIdμ to get astrophysical flux
    #this is likely not the most efficient way to do this.
    F = 2π * [Korg.trapezoid_rule(μ_surface_grid, μ_surface_grid .* I) for I in eachrow(I)]
    I, F
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
    trapezoid_rule(xs, fs)

Approximate the integral f(x) with the trapezoid rule over x-values `xs` given f(x) values `fs`.
"""
function trapezoid_rule(xs, fs)
    Δs = diff(xs)
    weights = [0 ; Δs] + [Δs ; 0]
    sum(0.5 * weights .* fs)
end

"""
    cumulative_trapezoid_rule(xs, fs)

Approximate the closed integral of f(x) from the first element of `xs` to each element of `xs` with
the trapezoid rule, given f values `fs`. Returns a vector of values of the same length as `xs` and 
`fs`.
"""
function cumulative_trapezoid_rule(xs, fs)
    [0 ; cumsum(0.5(fs[1:end-1] + fs[2:end]) .* diff(xs))]
end

"""
    exponential_integral_2(x)

Approximate second order exponential integral, E_2(x).  This stiches together several series 
expansions to get an approximation which is accurate within 1% for all `x`.
"""
function exponential_integral_2(x) 
    if x == 0
        0.0
    elseif x < 1.1
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
