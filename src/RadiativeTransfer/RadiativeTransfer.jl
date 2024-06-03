module RadiativeTransfer
using ...Korg: PlanarAtmosphere, ShellAtmosphere, CubicSplines, get_tau_5000s, get_zs

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

#TODO elliminate
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
function compute_astrophysical_flux(atm::PlanarAtmosphere, α, S, n_μ_points; include_inward_rays=false,
                                    τ_scheme="linear", I_scheme="linear_flux_only", α_ref=nothing)
    τ_ref = if !isnothing(α_ref) 
        get_tau_5000s(atm)
    else
        nothing
    end

    depths = get_zs(atm)
    depths .-= 2 * depths[end] #shift coordinates to avoid crossing zero

    radiative_transfer(α, S, depths, n_μ_points, include_inward_rays, false; 
                       α_ref=α_ref, τ_ref=τ_ref, I_scheme=I_scheme, τ_scheme=τ_scheme)
end
function compute_astrophysical_flux(atm::ShellAtmosphere, α, S, n_μ_points; include_inward_rays=false,
                                    τ_scheme=:anchored, I_scheme=:linear, α_ref=nothing)
    radii = [atm.R + l.z for l in atm.layers]
    photosphere_correction = radii[1]^2 / atm.R^2 

    τ_ref = if !isnothing(α_ref) #?
        get_tau_5000s(atm)
    else
        nothing
    end
    F, I = radiative_transfer(α, S, radii, n_μ_points, include_inward_rays, true; 
                              α_ref=α_ref, τ_ref=τ_ref, I_scheme=I_scheme, τ_scheme=τ_scheme)
    photosphere_correction .* F, I
end

"""
TODO

# Returns
- a spatial coordinate TODO increasing/decreasing along the ray
- the index of the lowest atmospheric layer pierced by each ray
"""
function calculate_rays(μ_surface_grid, spatial_coord, spherical)
    if spherical
        map(μ_surface_grid) do μ_surface
            b = spatial_coord[1] * sqrt(1 - μ_surface^2) # impact parameter of ray

            if b < spatial_coord[end] # ray goes below the atmosphere
                @. sqrt(spatial_coord^2 - b^2) 
            else
                # doing this with `findfirst` is messier at first and last index
                lowest_layer_index = argmin(abs.(spatial_coord .- b)) 
                if spatial_coord[lowest_layer_index] < b
                    lowest_layer_index -= 1
                end
                @. sqrt(spatial_coord[1:lowest_layer_index]^2 - b^2)
            end
        end
    else
        map(μ_surface_grid) do μ_surface
            spatial_coord ./ μ_surface
        end
    end
end

function radiative_transfer(α, S, spatial_coord, n_μ_points, include_inward_rays, spherical;
                            α_ref=nothing, τ_ref=nothing, I_scheme="linear_flux_only", τ_scheme="anchored")
    if τ_scheme == "spline"
        println("Warning: the spline τ scheme is not bug-free and may fail.")
    end

    μ_surface_grid, μ_weights = generate_mu_grid(n_μ_points) # TODO move this to the caller?

    # vector of path_length, layer_inds pairs
    rays = calculate_rays(μ_surface_grid, spatial_coord, spherical)

    # do inward rays either for everything, or just for the rays where we need to seed the bottom of
    # of the atmosphere
    inward_μ_surface_grid = if include_inward_rays
        -μ_surface_grid
    else
        -μ_surface_grid[length.(rays) .< length(spatial_coord)]
    end
    n_inward_rays = length(inward_μ_surface_grid)

    #type with which to preallocate arrays (enables autodiff)
    el_type = typeof(promote(spatial_coord[1], α[1], S[1], μ_surface_grid[1])[1])
    # intensity at every for every μ, λ, and layer. This is returned.
    # initialize with zeros because not every ray will pass through every layer
    I = if I_scheme == "linear_flux_only"
        # no "layers" dimension if we're only calculating the flux at the top of the atmosphere
        zeros(el_type, (n_inward_rays + length(μ_surface_grid), size(α, 2)))
    else
        zeros(el_type, (n_inward_rays + length(μ_surface_grid), size(α')...))
    end
    # preallocate a single τ vector which gets reused many times
    τ_buffer = Vector{el_type}(undef, length(spatial_coord)) 
    integrand_buffer = Vector{el_type}(undef, length(spatial_coord))
    log_τ_ref = log.(τ_ref) 
    #TODO precalculate λ-indenpendent quantities, at least for anchored τ, but maybe for other methods too

    # inward rays
    # TODO why is this twice as slow at the outward rays loop? (That performance is OK for now.)
    for μ_ind in 1:n_inward_rays
        path = -reverse(rays[μ_ind]) 
        layer_inds = length(path) : -1 : 1
        _radiative_transfer_core(μ_ind, layer_inds, n_inward_rays, path, τ_buffer, integrand_buffer, -log_τ_ref, α, S, I, spatial_coord, τ_ref, α_ref, τ_scheme, I_scheme)
    end

    # outward rays
    for μ_ind in n_inward_rays+1 : n_inward_rays+length(μ_surface_grid)
        path = rays[μ_ind - n_inward_rays]
        layer_inds = 1:length(path)
        _radiative_transfer_core(μ_ind, layer_inds, n_inward_rays, path, τ_buffer, integrand_buffer, log_τ_ref, α, S, I, spatial_coord, τ_ref, α_ref, τ_scheme, I_scheme)
    end

    #just the outward rays at the top layer
    surface_I = I[n_inward_rays+1:end, :, 1]
    F = 2π * (surface_I' * (μ_weights .* μ_surface_grid))

    F, I
end

function _radiative_transfer_core(μ_ind, layer_inds, n_inward_rays, path, τ_buffer, integrand_buffer, log_τ_ref, α, S, I, spatial_coord, τ_ref, α_ref, τ_scheme, I_scheme)
    # view into τ corresponding to the current ray TODO eliminate?
    τ = view(τ_buffer, layer_inds)

    integrand_factor = @. (spatial_coord[layer_inds] * τ_ref[layer_inds]) / (abs(path) * α_ref[layer_inds])

    for λ_ind in 1:size(α, 2)
        if length(path) <= 2
            # TODO try to do something smarter here
            # don't need to write anything because I is initialized to 0
            continue
        end

        # using more views below was not faster when I tested it
        # TODO: α is access in a cache-unfriendly way here
        if τ_scheme == "anchored"
            compute_tau_anchored!(τ, view(α, layer_inds, λ_ind), integrand_factor, log_τ_ref[layer_inds], integrand_buffer)
        elseif τ_scheme == "bezier"
            compute_tau_bezier!(τ, path, view(α, layer_inds, λ_ind))
        elseif τ_scheme == "spline"
            @info "spline scheme sometimes fails" #TODO 
            compute_tau_spline_analytic!(τ, path, view(α, layer_inds, λ_ind))
        else
            throw(ArgumentError("τ_scheme must be one of \"anchored\", \"bezier\", or \"spline\" (not recommended)"))
        end
        
        # these views into I are required because the function modifies I in place
        if I_scheme == "linear"
            #TODO switch S index order?
            linear_ray_transfer_integral!(view(I, μ_ind, λ_ind,  layer_inds), τ,
                                          view(S, layer_inds, λ_ind))
        elseif I_scheme == "linear_flux_only"
            # += because the intensity at the bottom of the atmosphere is already set for some rays
            I[μ_ind, λ_ind] += linear_ray_transfer_integral_flux_only(τ, view(S, layer_inds, λ_ind))
        elseif I_scheme == "bezier"
            bezier_ray_transfer_integral!(view(I, μ_ind, λ_ind, layer_inds), τ,
                                          view(S, layer_inds, λ_ind))
        else
            throw(ArgumentError("I_scheme must be one of \"linear\", \"bezier\", or \"linear_flux_only\""))
        end

        # set the intensity of the corresponding outward ray at the bottom of the atmosphere
        # this isn't correct for rays which go below the atmosphere, but the effect is immeasurable
        if μ_ind <= n_inward_rays # if ray is inwards
            if I_scheme == "linear_flux_only"
                # exp(-τ_buffer[1]) is the optical depth of the bottom of the atmosphere/end of the ray
                I[μ_ind + n_inward_rays, λ_ind] = I[μ_ind, λ_ind] * exp(-τ_buffer[1])
            else
                # TODO this branch does effectively nothing.  Which one is right?
                I[μ_ind + n_inward_rays, λ_ind, length(path)] = I[μ_ind, λ_ind, length(path)] 
            end
        end
    end
end

function compute_tau_anchored!(τ, α, integrand_factor, log_τ_ref, integrand_buffer)
    for k in eachindex(integrand_factor) #I can't figure out how to write this as a fast one-liner
        integrand_buffer[k] = α[k] * integrand_factor[k]
    end
    #MoogStyleTransfer.cumulative_trapezoid_rule!(τ, log_τ_ref, integrand_buffer) #TODO
    τ[1] = 0.0
    for i in 2:length(log_τ_ref)
        τ[i] = τ[i-1] + 0.5*(integrand_buffer[i]+integrand_buffer[i-1])*(log_τ_ref[i]-log_τ_ref[i-1])
    end
end

"""
    compute_tau_bezier(τ, s, α)

Compute optical depth (write to τ) along a ray with coordinate s and absorption coefficient α.  This 
is the method proposed in 
[de la Cruz Rodríguez and Piskunov 2013](https://ui.adsabs.harvard.edu/abs/2013ApJ...764...33D/abstract),
but the 
"""
function compute_tau_bezier!(τ, s, α)
    @assert length(τ) == length(s) == length(α) # because of the @inbounds below
    # how to get non-0 tau at first layer?
    τ[1] = 1e-5 
    C = fritsch_butland_C(s, α)
    # needed for numerical stability.  Threre is likely a smarter way to do this.
    clamp!(C, 1/2 * minimum(α), 2 * maximum(α))
    for i in 2:length(α)
        @inbounds τ[i] = τ[i-1] + (s[i-1] - s[i])/3 * (α[i] + α[i-1] + C[i-1])
    end
    ;
end
 
"""
    compute_tau_spline_analytic(τ, s, α)

Compute τ using a cubic spline to interpolate α.  This is not used by [`radiative_transfer`](@ref),
but is included for completeness.
"""
function compute_tau_spline_analytic!(τ, s, α)
    s = -s
    α_itp = CubicSplines.CubicSpline(s, α; extrapolate=true)
    CubicSplines.cumulative_integral!(τ, α_itp, s[1], s[end])
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
    @assert length(I) == length(τ) == length(S) # because of the @inbounds below

    I[end] = 0
    if length(τ) == 1
        return
    end

    for k in length(τ)-1:-1:1
        @inbounds δ = τ[k+1] - τ[k]
        @inbounds m = (S[k+1] - S[k])/δ
        @inbounds I[k] = (I[k+1] - S[k] -  m*(δ+1)) * (@fastmath exp(-δ)) + m + S[k]
    end
    ;
end

"""
TODO
"""
function linear_ray_transfer_integral_flux_only(τ, S)
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
    ray_transfer_integral!(I, τ, S)

TODO

Given τ and S along a ray (at a particular wavelength), compute the intensity at the end of the ray 
(the surface of the star).  This uses the method from 
[de la Cruz Rodríguez and Piskunov 2013](https://ui.adsabs.harvard.edu/abs/2013ApJ...764...33D/abstract).
"""
function bezier_ray_transfer_integral!(I, τ, S)
    @assert length(I) == length(τ) == length(S) # because of the @inbounds below
    I[end] = 0
    if length(τ) <= 1 
        return
    end

    C = fritsch_butland_C(τ, S)
    for k in length(τ)-1:-1:1
        @inbounds δ = τ[k+1] - τ[k]
        α = (2 + δ^2 - 2*δ - 2*exp(-δ)) / δ^2
        β = (2 - (2 + 2δ + δ^2)*exp(-δ)) / δ^2
        γ = (2*δ - 4 + (2δ + 4)*exp(-δ)) / δ^2

        @inbounds I[k] = I[k+1]*exp(-δ) + α*S[k] + β*S[k+1] + γ*C[k]
    end
    @inbounds I[1] *= exp(-τ[1]) #the second term isn't in the paper but it's necessary if τ[1] != 0
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