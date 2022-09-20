using FastGaussQuadrature: gausslegendre

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
function radiative_transfer(atm::PlanarAtmosphere, α, S, α_ref, mu_grid=nothing; bezier=false)
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

#TODO don't take tau ref here and in spherical transfer
function planar_transfer(α, S, z, τ_ref)
end

"""
    spherical_transfer(α, S, τ_ref, α_ref, radii, μ_surface_grid)

Perform radiative transfer along rays emerging at the μ values in `μ_surface_grid` in a spherically
symmetric atmosphere resolved at radii `radii` [cm]. See [`radiative_transfer`](@ref) for an 
explantion of the arguments. Note that `radii` should be in decreasing order.

Returns `(flux, intensity)`, where `flux` is the astrophysical flux, and `intensity`, a matrix of 
shape (mu values × wavelengths), is the surface intensity as a function of μ.
"""
function spherical_transfer(α, S, τ_ref, α_ref, radii, μ_surface_grid)
    μ_surface_grid, mu_weights = gausslegendre(length(μ_surface_grid))
    μ_surface_grid = @. μ_surface_grid/2 + 0.5
    mu_weights ./= 2

    #type with which to preallocate arrays (enables autodiff)
    el_type = typeof(promote(radii[1], α[1], S[1], μ_surface_grid[1])[1])

    # first calculate the wavelength-independent quantities:
    #  - ds/dr, the geometric path-length factor (n_layers × n_μ)
    #  - the index of the lowest atmospheric layer pierced by each ray
    ds_dr = Matrix{el_type}(undef, size(α, 1), length(μ_surface_grid))
    ds = Matrix{el_type}(undef, size(α, 1), length(μ_surface_grid))
    lowest_layer_indices = Vector{Int}(undef, length(μ_surface_grid))
    for (μ_ind, μ_surface) in enumerate(μ_surface_grid)
        b = radii[1] * sqrt(1 - μ_surface^2) # impact parameter of ray
        
        #doing this with `findfirst` is messier at first and last index
        i = argmin(abs.(radii .- b)) 
        if radii[i] < b
            i -= 1
        end
        lowest_layer_indices[μ_ind] = i

        ds_dr[1:i, μ_ind] = @. radii[1:i] ./ sqrt(radii[1:i]^2 - b^2) 
        #TODO try calculating ds directly
    end

    # iterate over λ in the outer loop, μ in the inner loop, calculating τ(r), then I(surface)
    I = Matrix{el_type}(undef, length(μ_surface_grid), size(α, 2)) #surface intensity (n_μ × n_λ)
    α_ds_dr = Vector{el_type}(undef, size(α, 1))  #integrand of τ integral, α * ds/dr
    τ_λ = Vector{el_type}(undef, size(α, 1)) 
    for λ_ind in 1:size(α, 2), μ_ind in 1:length(μ_surface_grid) 
        i = lowest_layer_indices[μ_ind]
        if i <= 2 #TODO ?
            I[μ_ind, λ_ind] = 0
            continue
        end
        for k in 1:i #I can't figure out how to write this as a fast one-liner
            α_ds_dr[k] = α[k, λ_ind] * ds_dr[k, μ_ind]
        end
        compute_tau_bezier!(view(τ_λ, 1:i), view(radii, 1:i), view(α_ds_dr,1:i))
        #if λ_ind == 1 && μ_ind == length(μ_surface_grid)
        #    display([τ_λ τ_ref])
        #end
        I[μ_ind, λ_ind] = ray_transfer_integral(view(τ_λ, 1:i), view(S, 1:i, λ_ind))

        # What happens at the lower boundery. Do we integrate through to the back of the star or 
        # stop? This branch could be factored out of this loop, which might speed things up.
        if i < length(radii)
            #if the ray never leaves the model atmosphere, include the contribution from the 
            #other side of the star.  Calculate the optical depths you get by reversing the τ_λ.
            #This is less accurate than actually integrating to find τ, but the effect is small.

            #This should probably by audited for off-by-one errors.  It's also inneficient.
            #TODO USE BEZIER SCHEME!
            #τ_prime = τ_λ[i] .+ cumsum(reverse(diff(view(τ_λ,1:i))))
            #I[μ_ind, λ_ind] += ray_transfer_integral(τ_prime, view(S,1:i-1,λ_ind))
        else 
            # otherwise assume I=S at atmosphere lower boundary.  This is a _tiny_ effect.
            I[μ_ind, λ_ind] += exp(-τ_λ[end]) * S[end, λ_ind]
        end
    end
    #calculate 2π∫μIdμ to get astrophysical flux
    #TODO - do this better.  Using a 5x finer mu grid results in a few 0.01% change in flux
    #F = 2π * [Korg.trapezoid_rule(μ_surface_grid, μ_surface_grid .* I) for I in eachcol(I)]
    println(size(mu_weights))
    println(size(I))
    println(size(I' * mu_weights))
    F = 2π * (I' * (mu_weights .* μ_surface_grid))
    I, F
end

#alpha should be increasing, s decreasing
#TODO off-by-one error, how to get non-0 tau at first layer
"""
TODO
"""
function compute_tau_bezier!(τ, s, α)
    τ[1] = 0
    C = fritsch_butland_C(s, α)
    for i in 2:length(α)
        τ[i] = τ[i-1] + (s[i-1] - s[i])/3 * (α[i] + α[i-1] + C[i-1])
    end
    ;
end

"""
TODO
"""
function ray_transfer_integral(τ, S)
    @assert length(τ) == length(S)
    if length(τ) == 1
        return 0.0
    end
    I = 0.0

    C = fritsch_butland_C(τ, S)
    for k in length(τ)-1:-1:1
        δ = τ[k+1] - τ[k]
        α = (2 + δ^2 - 2*δ - 2*exp(-δ)) / δ^2
        β = (2 - (2 + 2δ + δ^2)*exp(-δ)) / δ^2
        γ = (2*δ - 4 + (2δ + 4)*exp(-δ)) / δ^2
    
        I = I*exp(-δ) + α*S[k] + β*S[k+1] + γ*C[k]
    end
    I
end

"""
TODO
"""
function fritsch_butland_C(x, y)
    h = diff(x) #h[k] = x[k+1] - x[k]
    α = @. 1/3 * (1 + h[2:end]/(h[2:end] + h[1:end-1])) #α[k] is wrt h[k] and h[k-1]
    d = @. (y[2:end] - y[1:end-1])/h #d[k] is dₖ₊₀.₅ in paper
    yprime = @. (d[1:end-1] * d[2:end]) / (α*d[2:end] + (1-α)*d[1:end-1])

    C0 = @. y[2:end-1] + h[1:end-1]*yprime/2
    C1 = @. y[2:end-1] - h[2:end]*yprime/2

    ([C0 ; C1[end]] .+ [C0[1] ; C1]) ./ 2
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

