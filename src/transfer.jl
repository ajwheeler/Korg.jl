using FastGaussQuadrature: gausslegendre

"""
    radiative_transfer(atm::ModelAtmosphere, α, S, n_μ_points)

Returns the astrophysical flux at each wavelength.

inputs:
- `atm`: the model atmosphere.
- `α`: a matrix (atmospheric layers × wavelengths) containing the absoprtion coefficient
- `S`: the source fuction as a matrix of the same shape.
   rescale the total absorption to match the model atmosphere. This value should be calculated by 
   Korg.
- `n_μ_points`: the number of quadrature points to use when integrating over I_surface(μ) to obtain 
   the astrophysical flux.
"""
function radiative_transfer(atm::PlanarAtmosphere, α, S, n_μ_points; α5=nothing)
    zs = [l.z for l in atm.layers]
    planar_transfer(α, S, zs, n_μ_points, [l.tau_5000 for l in atm.layers], α5)
end
function radiative_transfer(atm::ShellAtmosphere, α, S, n_μ_points)
    radii = [atm.R + l.z for l in atm.layers]
    photosphere_correction = radii[1]^2 / atm.R^2 
    F, I = spherical_transfer(α, S, radii, n_μ_points)
    photosphere_correction .* F, I
end

"""
TODO
"""
function planar_transfer(α, S, z, n_μ_points, τ5, α5)
    μ_grid, μ_weights = generate_mu_grid(n_μ_points)

    #type with which to preallocate arrays (enables autodiff)
    el_type = typeof(promote(z[1], α[1], S[1], μ_grid[1])[1])

    #calculate adjustment
    korg_τ5 = Vector{el_type}(undef, size(α, 1))
    compute_tau_bezier!(korg_τ5, z, α5)
    spline_τ5 = Vector{el_type}(undef, size(α, 1))
    compute_tau_spline!(spline_τ5, z, α5)
    analytic_spline_τ5 = Vector{el_type}(undef, size(α, 1))
    compute_tau_spline_analytic!(analytic_spline_τ5, z, α5)
    #correction = τ5 ./ korg_τ5
    #display(correction)
    display([τ5 korg_τ5 spline_τ5 analytic_spline_τ5])

    # iterate over λ in the outer loop, μ in the inner loop, calculating τ(r), then I(surface)
    I = Matrix{el_type}(undef, length(μ_grid), size(α, 2)) #surface intensity (n_μ × n_λ)
    τ_λ = Vector{el_type}(undef, size(α, 1)) 
    for λ_ind in 1:size(α, 2), μ_ind in eachindex(μ_grid) 
        μ = μ_grid[μ_ind]

        compute_tau_bezier!(τ_λ, z ./ μ, view(α, :, λ_ind))
        #compute_tau_spline_analytic!(τ_λ, z ./ μ, view(α, :, λ_ind))
        #τ_λ .*= correction
        #τ_λ = τ5 ./ μ
        I[μ_ind, λ_ind] = ray_transfer_integral(τ_λ, view(S, :, λ_ind))

        # assume I=S at atmosphere lower boundary.  This is a _tiny_ effect.
        I[μ_ind, λ_ind] += exp(-τ_λ[end]) * S[end, λ_ind]
    end

    #calculate 2π∫μIdμ to get astrophysical flux
    F = 2π * (I' * (μ_weights .* μ_grid)) 
    F, I
end

"""
    spherical_transfer(α, S, τ_ref, α_ref, radii, μ_surface_grid)

Perform radiative transfer along rays emerging at the μ values in `μ_surface_grid` in a spherically
symmetric atmosphere resolved at radii `radii` [cm]. See [`radiative_transfer`](@ref) for an 
explantion of the arguments. Note that `radii` should be in decreasing order.

Returns `(flux, intensity)`, where `flux` is the astrophysical flux, and `intensity`, a matrix of 
shape (mu values × wavelengths), is the surface intensity as a function of μ.
"""
function spherical_transfer(α, S, radii, n_μ_points)
    μ_surface_grid, μ_weights = generate_mu_grid(n_μ_points)

    #type with which to preallocate arrays (enables autodiff)
    el_type = typeof(promote(radii[1], α[1], S[1], μ_surface_grid[1])[1])

    # first calculate the wavelength-independent quantities:
    #  - l, the coordinate along the ray
    #  - the index of the lowest atmospheric layer pierced by each ray
    l = Matrix{el_type}(undef, size(α, 1), length(μ_surface_grid))
    lowest_layer_indices = Vector{Int}(undef, length(μ_surface_grid))
    for (μ_ind, μ_surface) in enumerate(μ_surface_grid)
        b = radii[1] * sqrt(1 - μ_surface^2) # impact parameter of ray
        
        #doing this with `findfirst` is messier at first and last index
        i = argmin(abs.(radii .- b)) 
        if radii[i] < b
            i -= 1
        end

        lowest_layer_indices[μ_ind] = i
        l[1:i, μ_ind] = @. sqrt(radii[1:i]^2 - b^2)
    end

    # iterate over λ in the outer loop, μ in the inner loop, calculating τ(r), then I(surface)
    I = Matrix{el_type}(undef, length(μ_surface_grid), size(α, 2)) #surface intensity (n_μ × n_λ)
    τ_λ = Vector{el_type}(undef, size(α, 1)) 
    for λ_ind in 1:size(α, 2), μ_ind in 1:length(μ_surface_grid) 
        i = lowest_layer_indices[μ_ind]
        if i <= 2 #TODO ?
            I[μ_ind, λ_ind] = 0
            continue
        end
        compute_tau_bezier!(view(τ_λ, 1:i), view(l, 1:i, μ_ind), view(α, 1:i, λ_ind))
        I[μ_ind, λ_ind] = ray_transfer_integral(view(τ_λ, 1:i), view(S, 1:i, λ_ind))

        # At the lower boundary, we either integrate through to the back of the star or stop. 
        # This could be factored out of this loop, which might speed things up.
        if i < length(radii)
            # if the ray never leaves the model atmosphere, include the contribution from the 
            # other side of the star.
            
            # could preallocate for efficiency (make τ_λ one bigger to hold reversed tau)
            l_prime = [l[i, μ_ind] ; -view(l, i:-1:1, μ_ind)]
            α_prime = [α[i, λ_ind] ; view(α, i:-1:1, λ_ind)]
            τ_prime = similar(α_prime)
            compute_tau_bezier!(τ_prime, l_prime, α_prime)
            τ_prime .+= τ_λ[i]
            S_prime = [S[i, λ_ind] ; view(S,i:-1:1,λ_ind)]
            I[μ_ind, λ_ind] += ray_transfer_integral(τ_prime, S_prime)
        else 
            # otherwise assume I=S at atmosphere lower boundary.  This is a _tiny_ effect.
            I[μ_ind, λ_ind] += exp(-τ_λ[end]) * S[end, λ_ind]
        end
    end

    #calculate 2π∫μIdμ to get astrophysical flux
    F = 2π * (I' * (μ_weights .* μ_surface_grid))
    F, I
end

"""
TODO
"""
function generate_mu_grid(n_points)
    μ_grid, μ_weights = gausslegendre(n_points)
    μ_grid = @. μ_grid/2 + 0.5
    μ_weights ./= 2
    μ_grid, μ_weights
end

"""
TODO
"""
function compute_tau_bezier!(τ, s, α)
    #TODO off-by-one error, how to get non-0 tau at first layer
    τ[1] = 1e-5 #TODO
    C = fritsch_butland_C(s, α)
    clamp!(C, 1/2 * minimum(α), 2 * maximum(α)) # TODO needed for numerical stability?
    for i in 2:length(α)
        τ[i] = τ[i-1] + (s[i-1] - s[i])/3 * (α[i] + α[i-1] + C[i-1])
    end
    ;
end

function compute_tau_spline_analytic!(τ, s, α)
    s = -s
    α_itp = CubicSplines.CubicSpline(s, α; extrapolate=true)
    τ .= map(s) do layer_s
        CubicSplines.integral(α_itp, s[1], layer_s)
    end
end

"""
TODO
"""
function ray_transfer_integral(τ, S)
    @assert length(τ) == length(S)
    if length(τ) <= 1 #TODO?
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
    I * exp(-τ[1]) #the second term isn't in the paper but it's necessary if τ[1] != 0
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