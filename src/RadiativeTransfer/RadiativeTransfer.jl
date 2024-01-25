module RadiativeTransfer
using ...Korg: PlanarAtmosphere, ShellAtmosphere

# for generate_mu_grid
using FastGaussQuadrature: gausslegendre
"""
    generate_mu_grid(n_points)

Used by both radiative transfer schemes to compute quadature over μ. Returns `(μ_grid, μ_weights)`.
"""
function generate_mu_grid(n_points)
    μ_grid, μ_weights = gausslegendre(n_points)
    μ_grid = @. μ_grid/2 + 0.5
    μ_weights ./= 2
    μ_grid, μ_weights
end

include("BezierTransfer.jl")
include("MoogStyleTransfer.jl")

"""

# Arguments:
- `atm`: the model atmosphere.
- `α`: a matrix (atmospheric layers × wavelengths) containing the absoprtion coefficient
- `S`: the source fuction as a matrix of the same shape.
   rescale the total absorption to match the model atmosphere. This value should be calculated by 
   Korg.
- `μs`: the number of quadrature points to use when integrating over I_surface(μ) to obtain 
   the astrophysical flux. (TODO make this either a number or a vector of μ values)
"""
function compute_astrophysical_flux(atm::PlanarAtmosphere, α, S, n_μ_points; 
                                    tau_method=:anchored, I_method=:linear)
end
function compute_astrophysical_flux(atm::ShellAtmosphere, α, S, n_μ_points; do_negative_rays=false,
                                    tau_method=:anchored, I_method=:linear)
    radii = [atm.R + l.z for l in atm.layers]
    photosphere_correction = radii[1]^2 / atm.R^2 

    F, I = spherical_transfer(α, S, radii, n_μ_points, do_negative_rays)
    photosphere_correction .* F, I
end

"""
TODO

# Returns
- l, the coordinate along the ray
 - the index of the lowest atmospheric layer pierced by each ray
"""
function calculate_rays(μ_surface_grid, radii)
    map(enumerate(μ_surface_grid)) do (μ_ind, μ_surface)
        b = radii[1] * sqrt(1 - μ_surface^2) # impact parameter of ray
        

        if b < radii[end] # ray goes below the atmosphere
            @. sqrt(radii^2 - b^2)
        else
            #doing this with `findfirst` is messier at first and last index
            lowest_layer_index = argmin(abs.(radii .- b)) 
            if radii[lowest_layer_index] < b
                lowest_layer_index -= 1
            end
            path = @. sqrt(radii[1:lowest_layer_index]^2 - b^2)
            path[end] *= 2
            path = [path ; -path[end-1:-1:1]]
            path[1: min(length(radii), length(path))]
        end
    end
end

function spherical_transfer(α, S, radii, n_μ_points, do_negative_rays;
                            τ_ref=nothing, α_ref=nothing)
    μ_surface_grid, μ_weights = generate_mu_grid(n_μ_points)

    path_lengths = calculate_rays(μ_surface_grid, radii)

    # all_μ_surface_grid is the μ_surface_grid, but with negative μ values appended if they are 
    # being used
    all_μ_surface_grid = if do_negative_rays
        [μ_surface_grid ; -μ_surface_grid]
    else
        μ_surface_grid
    end

    #type with which to preallocate arrays (enables autodiff)
    el_type = typeof(promote(radii[1], α[1], S[1], μ_surface_grid[1])[1])

    #TODO precalculate λ-indenpendent quantities, at least for anchored τ, but maybe for other methods too

    # intensity at every layer, for every μ, for every λ. This is returned
    I = Array{el_type}(undef, length(all_μ_surface_grid), size(α)...) 
    # preallocate a single τ vector which gets reused many times
    τ_λ = Vector{el_type}(undef, length(radii)) 
    for μ_ind in 1:length(all_μ_surface_grid)
        # index of the ray in the μ_surface_grid with the same |μ| as this ray
        positive_μ_ind = μ_ind <= length(μ_surface_grid) ? μ_ind : μ_ind - length(μ_surface_grid)

        path = path_lengths[positive_μ_ind]
        lowest_layer_ind = length(path)

        # indices of layers along this ray, in the correct order
        layer_inds = if μ_ind <= length(μ_surface_grid)
            1:lowest_layer_ind #ray coming out
        else
            path .*= -1
            lowest_layer_ind:-1:1 #ray going in
        end

        for λ_ind in 1:size(α, 2)
            if lowest_layer_ind <= 2
                I[μ_ind, :, λ_ind] .= 0
                continue
            end

            # TODO switch this to whatever
            BezierTransfer.compute_tau_bezier!(view(τ_λ, layer_inds),
                                               path[layer_inds],
                                               view(α, layer_inds, λ_ind))

            # TODO switch this to whatever
            bezier_ray_transfer_integral!(view(I, μ_ind, layer_inds, λ_ind),
                                   view(τ_λ, layer_inds),
                                   view(S, layer_inds, λ_ind))

            if λ_ind == 4000 && μ_ind in [3, 23]
                println("μ_ind = $μ_ind, μ = $(all_μ_surface_grid[μ_ind])")
                display(["path" "alpha" "τ" "I" "S" ; 
                         path view(α, layer_inds, λ_ind) view(τ_λ, layer_inds) view(I, μ_ind, layer_inds, λ_ind)  view(S, layer_inds, λ_ind)])
                println()
                println()
            end
        end
    end

    # TODO how to account for this nicely?
    #for μ_ind in eachindex(μ_surface_grid)
    #    i = lowest_layer_indices[μ_ind]
    #    for λ_ind in 1:size(α, 2)
    #        if (i < length(radii)) && (τ_λ[i] < 100.0)
    #            #if the ray never leaves the model atmosphere, include the contribution from the 
    #            #other side of the star.  We also check that τ is not ridiculously large, because
    #            #numerical precision can lead to adjacent layers getting equal τ on the other side of 
    #            #the star
    #            
    #            # could preallocate for efficiency (make τ_λ one bigger to hold reversed tau)
    #            l_prime = [l[i, μ_ind] ; -view(l, i:-1:1, μ_ind)]
    #            α_prime = [α[i, λ_ind] ; view(α, i:-1:1, λ_ind)]
    #            τ_prime = similar(α_prime)

    #            #TODO switch this to whatever
    #            compute_tau_bezier!(τ_prime, l_prime, α_prime)
    #            #compute_tau_spline_analytic!(τ_prime, l_prime, α_prime)
    #            τ_prime .+= τ_λ[i]
    #            S_prime = [S[i, λ_ind] ; view(S,i:-1:1,λ_ind)]

    #            # TODO switch this to whatever
    #            # TODO make this work
    #            #I[μ_ind, λ_ind] += ray_transfer_integral(τ_prime, S_prime)
    #        else 
    #            # otherwise assume I=S at atmosphere lower boundary.  This is a _tiny_ effect.
    #            I[μ_ind, λ_ind] += exp(-τ_λ[i]) * S[end, λ_ind]
    #        end
    #    end
    #end
    #calculate 2π∫μIdμ to get astrophysical flux

    #just the outward rays at the top layer
    surface_I = I[1:length(μ_surface_grid), 1, :]
    F = 2π * (surface_I' * (μ_weights .* μ_surface_grid))

    F, I
end

"""
    ray_transfer_integral!(I, τ, S)

TODO

Given τ and S along a ray (at a particular wavelength), compute the intensity at the end of the ray 
(the surface of the star).  This uses the method from 
[de la Cruz Rodríguez and Piskunov 2013](https://ui.adsabs.harvard.edu/abs/2013ApJ...764...33D/abstract).
"""
function bezier_ray_transfer_integral!(I, τ, S)
    @assert length(τ) == length(S)  == length(I)
    I[end] = 0
    if length(τ) <= 1 
        return
    end

    C = fritsch_butland_C(τ, S)
    for k in length(τ)-1:-1:1
        δ = τ[k+1] - τ[k]
        α = (2 + δ^2 - 2*δ - 2*exp(-δ)) / δ^2
        β = (2 - (2 + 2δ + δ^2)*exp(-δ)) / δ^2
        γ = (2*δ - 4 + (2δ + 4)*exp(-δ)) / δ^2
    
        I[k] = I[k+1]*exp(-δ) + α*S[k] + β*S[k+1] + γ*C[k]
    end
    I[1] *= exp(-τ[1]) #the second term isn't in the paper but it's necessary if τ[1] != 0
    ;
end

"""
    ray_transfer_integral(I, τ, S)

TODO

Compute exactly the solution to the transfer integral obtained be linearly interpolating the source 
function, `S` across optical depths `τ`, without approximating the factor of exp(-τ).

This breaks the integral into the sum of integrals of the form 
\$\\int (m\\tau + b) \\exp(-\\tau)\$ d\\tau\$ , 
which is equal to
\$ -\\exp(-\\tau) (m*\\tau + b + m)\$.
"""
function linear_ray_transfer_integral!(I, τ, S)
    @assert length(τ) == length(S)  == length(I)
    I[end] = 0
    if length(τ) == 1
        return
    end

    for k in length(τ)-1:-1:1
        δ = τ[k+1] - τ[k]
        m = (S[k+1] - S[k])/δ
        I[k] = (I[k+1] - S[k] -  m*(δ+1)) * exp(-δ) + m + S[k]
    end
    ;
end

"""
    fritsch_butland_C(x, y)

Given a set of x and y values, compute the bezier control points using the method of 
[Fritch & Butland 1984](https://doi.org/10.1137/0905021), as suggested in 
[de la Cruz Rodríguez and Piskunov 2013](https://ui.adsabs.harvard.edu/abs/2013ApJ...764...33D/abstract).
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

end #module