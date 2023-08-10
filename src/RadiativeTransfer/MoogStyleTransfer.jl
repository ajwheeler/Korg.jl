"""
Korg's default radiative transfer implementation. 
See also: `RadiativeTransfer.BezierTransfer`
"""
module MoogStyleTransfer 
using ..RadiativeTransfer: generate_mu_grid
using ...Korg: PlanarAtmosphere, ShellAtmosphere
    
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

Keyword arguments:
- `n_nu_points`: the number of rays to use when integrating over I_surface(μ) to obtain the 
   astrophysical flux.  This doesn't do anything if `atm` is a [`PlanarAtmosphere`](@ref).
"""
function radiative_transfer(atm::PlanarAtmosphere, α, S, α_ref, n_mu_points=nothing)
    τ5 = [l.tau_5000 for l in atm.layers] #τ at 5000 Å according to model atmosphere
    planar_transfer(α, S, τ5, α_ref), nothing
end
function radiative_transfer(atm::ShellAtmosphere, α, S, α_ref, n_mu_points)
    τ5 = [l.tau_5000 for l in atm.layers] #τ at 5000 Å according to model atmosphere
    radii = [atm.R + l.z for l in atm.layers]
    photosphere_correction = radii[1]^2 / atm.R^2
    #discard I, take F only
    F, I = spherical_transfer(α, S, τ5, α_ref, radii, n_mu_points)
    photosphere_correction * F, I
end

"""
    planar_transfer(α, S, τ_ref, α_ref)

Returns the astrophysical flux. See [`radiative_transfer`](@ref) for an explantion of the arguments.
"""
function planar_transfer(α, S, τ_ref, α_ref)
    log_τ_ref = log.(τ_ref)
    τ_ref_α_ref = τ_ref ./ α_ref
    I = Vector{eltype(α)}(undef, size(α, 2))
    τ_λ = Vector{eltype(α)}(undef, size(α, 1))
    for λ_ind in 1:size(α, 2)
        @inbounds cumulative_trapezoid_rule!(τ_λ, log_τ_ref, view(α, :, λ_ind) .* τ_ref_α_ref)
        I[λ_ind] = 2π * all_mu_transfer_integral(τ_λ, view(S, :, λ_ind))
    end
    I
end

"""
    spherical_transfer(α, S, τ_ref, α_ref, radii, μ_surface_grid)

Perform radiative transfer along rays emerging at the μ values in `μ_surface_grid` in a spherically
symmetric atmosphere resolved at radii `radii` [cm]. See [`radiative_transfer`](@ref) for an 
explantion of the arguments. Note that `radii` should be in decreasing order.

Returns `(flux, intensity)`, where `flux` is the astrophysical flux, and `intensity`, a matrix of 
shape (mu values × wavelengths), is the surface intensity as a function of μ.
"""
function spherical_transfer(α, S, τ_ref, α_ref, radii, n_μ_points)
    μ_surface_grid, μ_weights = generate_mu_grid(n_μ_points)

    R, r0 = radii[1], radii[end] #lower bound of atmosphere

    #precompute for use in transfer integral
    log_τ_ref = log.(τ_ref) 
    τ_ref_α_ref = τ_ref ./ α_ref

    #general element type with which to preallocate arrays
    el_type = typeof(promote(radii[1], α[1], S[1], μ_surface_grid[1])[1])

    #preallocations
    #geometric path-length correction, and other wavelength-indenpendent integrand factors
    integrand_factor = Matrix{el_type}(undef, size(α, 1), length(μ_surface_grid))
    #the index of the lowest atmospheric layer pierced by each ray
    lowest_layer_indices = Vector{Int}(undef, length(μ_surface_grid))
    for (μ_ind, μ_surface) in enumerate(μ_surface_grid)
        # impact parameter of ray
        b = R * sqrt(1 - μ_surface^2)
        
        #doing this with `findfirst` is messier at first and last index
        i = argmin(abs.(radii .- b)) 
        if radii[i] < b
            i -= 1
        end
        lowest_layer_indices[μ_ind] = i

        integrand_factor[1:i, μ_ind] = @. radii[1:i] ./ sqrt(radii[1:i]^2 - b^2) * τ_ref_α_ref[1:i]
    end

    #preallocations
    I = Matrix{el_type}(undef, length(μ_surface_grid), size(α, 2)) #the surface intensity
    integrand = Vector{el_type}(undef, size(α, 1))                 #integrand of τ integral 
    τ_λ = Vector{el_type}(undef, size(α, 1))                       #optical depth at a particular λ
    #iterate over λ in the outer loop, μ in the inner loop
    for λ_ind in 1:size(α, 2), μ_ind in 1:length(μ_surface_grid) 
        i = lowest_layer_indices[μ_ind]
        for k in 1:i #I can't figure out how to write this as a fast one-liner
            integrand[k] = α[k, λ_ind] * integrand_factor[k, μ_ind]
        end
        cumulative_trapezoid_rule!(τ_λ, log_τ_ref, integrand, i) #compute τ_λ
        I[μ_ind, λ_ind] = ray_transfer_integral(view(τ_λ,1:i), view(S,1:i,λ_ind))

        #this branch could be factored out of this loop, which might speed things up.
        if (i < length(radii)) && (τ_λ[i] < 100.0)
            #if the ray never leaves the model atmosphere, include the contribution from the 
            #other side of the star.  We also check that τ is not ridiculously large, because
            #numerical precision can lead to adjacent layers getting equal τ on the other side of 
            #the star
            
            #This is less accurate than actually integrating to find τ, but the effect is small.
            #This should probably by audited for off-by-one errors.  It's also inneficient.
            τ_prime = τ_λ[i] .+ [cumsum(reverse(diff(view(τ_λ,1:i)))) ; 0]
            I[μ_ind, λ_ind] += ray_transfer_integral(view(τ_prime, 1:i), view(S,1:i,λ_ind))
        else #otherwise assume I=S at atmosphere lower boundary.  This is a _tiny_ effect.
            I[μ_ind, λ_ind] += exp(-τ_λ[i]) * S[end, λ_ind]
        end
    end
    #calculate 2π∫μIdμ to get astrophysical flux
    F = 2π * (I' * (μ_weights .* μ_surface_grid))
    F, I
end

"""
    all_mu_transfer_integral(τ, S)

Compute exactly the solution to the transfer integral obtained be linearly interpolating the source 
function, `S` across optical depths `τ`, without approximating the factor of E₂(τ).
"""
function all_mu_transfer_integral(τ, S)
    I = 0
    for i in 1:length(τ)-1
        @inbounds m = (S[i+1] - S[i])/(τ[i+1] - τ[i])
        @inbounds b = S[i] - m*τ[i]
        @inbounds I += (_plane_parallel_approximate_transfer_integral(τ[i+1], m, b) - 
                        _plane_parallel_approximate_transfer_integral(τ[i], m, b))
    end
    I
end

"""
    ray_transfer_integral(τ, S)

Compute exactly the solution to the transfer integral obtained be linearly interpolating the source 
function, `S` across optical depths `τ`, without approximating the factor of exp(-τ).

This breaks the integral into the sum of integrals of the form 
\$\\int (m\\tau + b) \\exp(-\\tau)\$ d\\tau\$ , 
which is equal to
\$ -\\exp(-\\tau) (m*\\tau + b + m)\$.
"""
function ray_transfer_integral(τ, S)
    if length(τ) == 1
        return 0.0
    end
    I = 0.0
    next_exp_negτ = exp(-τ[1])
    for i in 1:length(τ)-1
        @inbounds m = (S[i+1] - S[i])/(τ[i+1] - τ[i])
        cur_exp_negτ = next_exp_negτ
        @inbounds next_exp_negτ = exp(-τ[i+1])
        @inbounds I += (-next_exp_negτ * (S[i+1] + m) + cur_exp_negτ * (S[i] + m))
    end
    I
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
    cumulative_trapezoid_rule!(out, xs, fs, [len=length(xs)])

Approximate the closed integral of f(x) from the first element of `xs` to each element of `xs` 
(up to `len`) with the trapezoid rule, given f values `fs`. 

Assigns to the preallocated vector `out`.
"""
function cumulative_trapezoid_rule!(out, xs, fs, len=length(xs))
    out[1] = 0.0
    for i in 2:len
        out[i] = out[i-1] + 0.5*(fs[i]+fs[i-1])*(xs[i]-xs[i-1])
    end
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

end # module