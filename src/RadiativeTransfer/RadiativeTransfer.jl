module RadiativeTransfer
using ...Korg: PlanarAtmosphere, ShellAtmosphere, get_tau_5000s, get_zs
using FastGaussQuadrature: gausslegendre

"""
    generate_mu_grid(n_points)
    generate_mu_grid(μ_values)

Used by both radiative transfer schemes to compute quadrature over μ. Returns `(μ_grid, μ_weights)`.
If an integer is passed, generate a grid of Gauss-Legendre quadrature points of corresponding size.
Otherwise, use the provided μ values.
"""
function generate_mu_grid(n_points::Integer)
    μ_grid, μ_weights = gausslegendre(n_points)
    μ_grid = @. μ_grid / 2 + 0.5
    μ_weights ./= 2
    μ_grid, μ_weights
end
function generate_mu_grid(μ_grid::AbstractVector{<:Real})
    if length(μ_grid) == 1
        return μ_grid, [1.0]
    end

    if !issorted(μ_grid) || μ_grid[1] < 0 || μ_grid[end] > 1
        throw(ArgumentError("μ_grid must be sorted and bounded between 0 and 1"))
    end
    Δ = diff(μ_grid)
    μ_weights = 0.5 * [Δ[1]; (@. Δ[1:end-1] + Δ[2:end]); Δ[end]]
    μ_grid, μ_weights
end

"""
    radiative_transfer(atm, α, S, n_μ_points; kwargs...)
    radiative_transfer(α, S, spatial_coord, n_μ_points, spherical; kwargs...)

# Arguments (note that some of these are only required for one method):

  - `atm`: the model atmosphere.
  - `α`: a matrix (atmospheric layers × wavelengths) containing the absorption coefficient
  - `S`: the source function as a matrix of the same shape.
    rescale the total absorption to match the model atmosphere. This value should be calculated by
    Korg.
  - `μ_points`: the number of quadrature points to use when integrating over I_surface(μ) to obtain
    the astrophysical flux, or, alternatively, the specific μ values to use (in which case the
    integral is done with the trapezoid rule).
  - `spatial_coord`: a physical distance coordinate (radius for a spherical atmosphere, height above
    the photosphere for a plane-parallel atmosphere).
  - `spherical`: whether the atmosphere is spherical or plane-parallel.

# Keyword Arguments:

  - `include_inward_rays` (default: `false`): if true, light propagating into the star (negative μs)
    is included.  If false, only those which are needed to seed the intensity at the bottom of the
    atmosphere are included.
  - `τ_scheme` (default: "anchored"): how to compute the optical depth.  Options are "linear" and
    "bezier" (not recommended).
  - `I_scheme` (default: "linear_flux_only"): how to compute the intensity.  Options are "linear",
    "linear_flux_only", and "bezier".  "linear_flux_only" is the fastest, but does not return the
    intensity values anywhere except at the top of the atmosphere.  "linear" performs an equivalent
    calculation, but stores the intensity at every layer.  "bezier" is not recommended.
"""
function radiative_transfer(atm::PlanarAtmosphere, α, S, μ_points;
                            α_ref=nothing, τ_ref=isnothing(α_ref) ? nothing : get_tau_5000s(atm),
                            τ_scheme="anchored", I_scheme="linear_flux_only", kwargs...)
    radiative_transfer(α, S, get_zs(atm), μ_points, false; α_ref=α_ref, τ_ref=τ_ref,
                       I_scheme=I_scheme, τ_scheme=τ_scheme, kwargs...)
end
function radiative_transfer(atm::ShellAtmosphere, α, S, μ_points;
                            α_ref=nothing, τ_ref=isnothing(α_ref) ? nothing : get_tau_5000s(atm),
                            τ_scheme="anchored", I_scheme="linear", kwargs...)
    radii = [atm.R + l.z for l in atm.layers]
    photosphere_correction = radii[1]^2 / atm.R^2
    F, others... = radiative_transfer(α, S, radii, μ_points, true; α_ref=α_ref, τ_ref=τ_ref,
                                      I_scheme=I_scheme, τ_scheme=τ_scheme, kwargs...)
    photosphere_correction .* F, others...
end
function radiative_transfer(α, S, spatial_coord, μ_points, spherical;
                            include_inward_rays=false,
                            α_ref=nothing, τ_ref=nothing, I_scheme="linear_flux_only",
                            τ_scheme="anchored")
    if I_scheme == "linear_flux_only" && τ_scheme == "anchored" && !spherical
        I_scheme = "linear_flux_only_expint"
        # in this special case, we can use exponential integral tricks
        μ_surface_grid, μ_weights = [1], [1]
    else
        μ_surface_grid, μ_weights = generate_mu_grid(μ_points)
    end

    # distance along ray, and derivative wrt spatial coord
    rays = calculate_rays(μ_surface_grid, spatial_coord, spherical)

    # do inward rays either for everything, or just for the rays where we need to seed the bottom of
    # of the atmosphere
    inward_μ_surface_grid = if include_inward_rays
        -μ_surface_grid
    else
        -μ_surface_grid[length.(first.(rays)).<length(spatial_coord)]
    end
    n_inward_rays = length(inward_μ_surface_grid)

    #type with which to preallocate arrays (enables autodiff)
    el_type = typeof(promote(spatial_coord[1], α[1], S[1], μ_surface_grid[1])[1])
    # intensity at every for every μ, λ, and layer. This is returned.
    # initialize with zeros because not every ray will pass through every layer
    I = if startswith(I_scheme, "linear_flux_only") # may or may not end in _expint
        # no "layers" dimension if we're only calculating the flux at the top of the atmosphere
        zeros(el_type, (n_inward_rays + length(μ_surface_grid), size(α, 2)))
    else
        zeros(el_type, (n_inward_rays + length(μ_surface_grid), size(α')...))
    end
    # preallocate a single τ vector which gets reused many times
    τ_buffer = Vector{el_type}(undef, length(spatial_coord))
    integrand_buffer = Vector{el_type}(undef, length(spatial_coord))
    log_τ_ref = log.(τ_ref)

    # inward rays (this twice as slow at the outward rays loop, which would be good to improve)
    for μ_ind in 1:n_inward_rays
        path, dsdz = reverse.(rays[μ_ind])
        layer_inds = length(path):-1:1
        radiative_transfer_core(μ_ind, layer_inds, n_inward_rays, -path, dsdz,
                                τ_buffer, integrand_buffer, -log_τ_ref, α, S, I,
                                τ_ref, α_ref, τ_scheme, I_scheme)
    end

    # outward rays
    for μ_ind in n_inward_rays+1:n_inward_rays+length(μ_surface_grid)
        path, dsdz = rays[μ_ind-n_inward_rays]
        layer_inds = 1:length(path)
        radiative_transfer_core(μ_ind, layer_inds, n_inward_rays, path, dsdz,
                                τ_buffer, integrand_buffer, log_τ_ref, α, S, I,
                                τ_ref, α_ref, τ_scheme, I_scheme)
    end

    #just the outward rays at the top layer
    surface_I = I[n_inward_rays+1:end, :, 1]
    F = 2π * (surface_I' * (μ_weights .* μ_surface_grid))

    F, I, μ_surface_grid, μ_weights
end

"""
Compute the intensity along a single ray.

n.b. this function has an additional I_scheme ("linear_flux_only_expint") that radiative_transfer
will automatically switch to when appropriate.
"""
function radiative_transfer_core(μ_ind, layer_inds, n_inward_rays, path, dsdz, τ_buffer,
                                 integrand_buffer, log_τ_ref, α, S, I, τ_ref, α_ref,
                                 τ_scheme, I_scheme)
    if length(path) == 1 && ((I_scheme == "bezier") || (τ_scheme == "bezier"))
        # these schemes requires two layers minimum
        I[μ_ind, :, 1] .= 0.0
        return
    end

    # view into τ corresponding to the current ray
    τ = view(τ_buffer, layer_inds)

    # this is τref/αref * ds/dz
    integrand_factor = @. τ_ref[layer_inds] / α_ref[layer_inds] * dsdz

    for λ_ind in 1:size(α, 2)
        # using more views below was not faster when I tested it
        # α is access in a cache-unfriendly way here. Fixing that might speex things up.
        if τ_scheme == "anchored"
            compute_tau_anchored!(τ, view(α, layer_inds, λ_ind), integrand_factor,
                                  log_τ_ref[layer_inds], integrand_buffer)
        elseif τ_scheme == "bezier"
            compute_tau_bezier!(τ, path, view(α, layer_inds, λ_ind))
        else
            throw(ArgumentError("τ_scheme must be one of \"anchored\" or \"bezier\""))
        end

        # these views into I are required because the function modifies I in place
        if I_scheme == "linear"
            # switching the S index order might speed things up
            compute_I_linear!(view(I, μ_ind, λ_ind, layer_inds), τ,
                              view(S, layer_inds, λ_ind))
        elseif I_scheme == "linear_flux_only"
            # += because the intensity at the bottom of the atmosphere is already set for some rays
            I[μ_ind, λ_ind] += compute_I_linear_flux_only(τ, view(S, layer_inds, λ_ind))
        elseif I_scheme == "linear_flux_only_expint"
            I[μ_ind, λ_ind] += compute_F_flux_only_expint(τ, view(S, layer_inds, λ_ind))
        elseif I_scheme == "bezier"
            compute_I_bezier!(view(I, μ_ind, λ_ind, layer_inds), τ,
                              view(S, layer_inds, λ_ind))
        else
            throw(ArgumentError("I_scheme must be one of \"linear\", \"bezier\", or \"linear_flux_only\""))
        end

        # set the intensity of the corresponding outward ray at the bottom of the atmosphere
        # this isn't correct for rays which go below the atmosphere, but the effect is immeasurable
        if μ_ind <= n_inward_rays # if ray is inwards
            if startswith(I_scheme, "linear_flux_only") # may or may not end in _expint
                I[μ_ind+n_inward_rays, λ_ind] = I[μ_ind, λ_ind] * exp(-τ[end])
            else
                I[μ_ind+n_inward_rays, λ_ind, length(path)] = I[μ_ind, λ_ind, length(path)]
            end
        end
    end
end

"""
    calculate_rays(μ_surface_grid, spatial_coord, spherical)

# Arguments

  - `μ_surface_grid`: the μ values at the surface of the star corresponding to the rays along which
    the optical depth and intensity will be calculated.
  - `spatial_coord`: a physical distance coordinate.  This is radius for a spherical atmosphere, and
    height above the photosphere for a plane-parallel atmosphere.

# Returns

A vector of pairs `(s, ds/dz)`, where `s` is the distance along the ray and z is the model
atmosphere spatial coordinate.
"""
function calculate_rays(μ_surface_grid, spatial_coord, spherical)
    if spherical # spatial_coord is radius
        map(μ_surface_grid) do μ_surface
            b = spatial_coord[1] * sqrt(1 - μ_surface^2) # impact parameter of ray

            lowest_layer_index = if b < spatial_coord[end] # ray goes below the atmosphere
                length(spatial_coord)
            else
                # doing this with `findfirst` is messier at first and last index
                lowest_layer_index = argmin(abs.(spatial_coord .- b))
                if spatial_coord[lowest_layer_index] < b
                    lowest_layer_index -= 1
                end
                lowest_layer_index
            end
            s = @. sqrt(spatial_coord[1:lowest_layer_index]^2 - b^2)
            dsdr = @. spatial_coord[1:lowest_layer_index] ./ s
            s, dsdr
        end
    else # spatial_coord measured relative to whatever
        map(μ_surface_grid) do μ_surface
            (spatial_coord ./ μ_surface), ones(length(spatial_coord)) ./ μ_surface
        end
    end
end

function compute_tau_anchored!(τ, α, integrand_factor, log_τ_ref, integrand_buffer)
    for k in eachindex(integrand_factor) #I can't figure out how to write this as a fast one-liner
        integrand_buffer[k] = α[k] * integrand_factor[k]
    end
    τ[1] = 0.0
    for i in 2:length(log_τ_ref)
        τ[i] = τ[i-1] +
               0.5 * (integrand_buffer[i] + integrand_buffer[i-1]) * (log_τ_ref[i] - log_τ_ref[i-1])
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
    clamp!(C, 1 / 2 * minimum(α), 2 * maximum(α))
    for i in 2:length(α)
        @inbounds τ[i] = τ[i-1] + (s[i-1] - s[i]) / 3 * (α[i] + α[i-1] + C[i-1])
    end
end

"""
    compute_I_linear!(I, τ, S)

Compute exactly the solution to the transfer integral obtained be linearly interpolating the source
function, `S` across optical depths `τ`, without approximating the factor of exp(-τ).

This breaks the integral into the sum of integrals of the form
\$\\int (m\\tau + b) \\exp(-\\tau)\$ d\\tau\$ ,
which is equal to
\$ -\\exp(-\\tau) (m*\\tau + b + m)\$.
"""
function compute_I_linear!(I, τ, S)
    @assert length(I) == length(τ) == length(S) # because of the @inbounds below

    if length(τ) == 1
        return
    end

    for k in length(τ)-1:-1:1
        @inbounds δ = τ[k+1] - τ[k]
        @inbounds m = (S[k+1] - S[k]) / δ
        @inbounds I[k] = (I[k+1] - S[k] - m * (δ + 1)) * (@fastmath exp(-δ)) + m + S[k]
    end
end

"""
    compute_I_linear_flux_only(τ, S)

Returns the intensity at the end of the ray (the surface of the star) given τ and S along the ray.
Uses the same numerical method as [`compute_I_linear!`](@ref), but doesn't retain the intensity at
each layer.
"""
function compute_I_linear_flux_only(τ, S)
    if length(τ) == 1
        return 0.0
    end
    I = 0.0
    next_exp_negτ = exp(-τ[1])
    for i in 1:length(τ)-1
        @inbounds Δτ = τ[i+1] - τ[i]
        # fix the case where large τ causes numerically 0 Δτ
        Δτ += (Δτ == 0) # if it's 0, make it 1
        @inbounds m = (S[i+1] - S[i]) / Δτ

        cur_exp_negτ = next_exp_negτ
        @inbounds next_exp_negτ = exp(-τ[i+1])
        @inbounds I += (-next_exp_negτ * (S[i+1] + m) + cur_exp_negτ * (S[i] + m))
    end
    I
end

"""
    compute_I_bezier!(I, τ, S)

Given τ and S along a ray (at a particular wavelength), compute the intensity at the end of the ray
(the surface of the star).  This uses the method from
[de la Cruz Rodríguez and Piskunov 2013](https://ui.adsabs.harvard.edu/abs/2013ApJ...764...33D/abstract).
"""
function compute_I_bezier!(I, τ, S)
    @assert length(I) == length(τ) == length(S) # because of the @inbounds below
    I[end] = 0
    if length(τ) <= 1
        return
    end

    C = fritsch_butland_C(τ, S)
    for k in length(τ)-1:-1:1
        @inbounds δ = τ[k+1] - τ[k]
        α = (2 + δ^2 - 2 * δ - 2 * exp(-δ)) / δ^2
        β = (2 - (2 + 2δ + δ^2) * exp(-δ)) / δ^2
        γ = (2 * δ - 4 + (2δ + 4) * exp(-δ)) / δ^2

        @inbounds I[k] = I[k+1] * exp(-δ) + α * S[k] + β * S[k+1] + γ * C[k]
    end
    @inbounds I[1] *= exp(-τ[1]) #the second term isn't in the paper but it's necessary if τ[1] != 0
end

"""
    fritsch_butland_C(x, y)

Given a set of x and y values, compute the bezier control points using the method of
[Fritch & Butland 1984](https://doi.org/10.1137/0905021), as suggested in
[de la Cruz Rodríguez and Piskunov 2013](https://ui.adsabs.harvard.edu/abs/2013ApJ...764...33D/abstract).

Used in [`compute_I_bezier!`](@ref).
"""
function fritsch_butland_C(x, y)
    h = diff(x) #h[k] = x[k+1] - x[k]
    α = @. 1 / 3 * (1 + h[2:end] / (h[2:end] + h[1:end-1])) #α[k] is wrt h[k] and h[k-1]
    d = @. (y[2:end] - y[1:end-1]) / h #d[k] is dₖ₊₀.₅ in paper
    yprime = @. (d[1:end-1] * d[2:end]) / (α * d[2:end] + (1 - α) * d[1:end-1])

    C0 = @. y[2:end-1] + h[1:end-1] * yprime / 2
    C1 = @. y[2:end-1] - h[2:end] * yprime / 2

    ([C0; C1[end]] .+ [C0[1]; C1]) ./ 2
end

"""
    compute_F_flux_only_expint(τ, S)

Compute the astrophysical flux, F, by linearly interpolating the source
function, `S` across optical depths `τ`.  Handle the integral over μ analytically using E₂.
"""
function compute_F_flux_only_expint(τ, S)
    I = 0
    for i in 1:length(τ)-1
        @inbounds m = (S[i+1] - S[i]) / (τ[i+1] - τ[i])
        @inbounds b = S[i] - m * τ[i]
        @inbounds I += (expint_transfer_integral_core(τ[i+1], m, b) -
                        expint_transfer_integral_core(τ[i], m, b))
    end
    I
end

"""
    expint_transfer_integral_core(τ, m, b)

The exact solution to \$\\int (m\\tau + b) E_2(\\tau)\$ d\\tau\$.
The exponential integral function, expint, captures the integral over the disk of the star to
get the emergent astrophysical flux. You can verify it by substituting the variable of integration
in the exponential integal, t, with mu=1/t.
"""
function expint_transfer_integral_core(τ, m, b)
    1 / 6 * (τ * exponential_integral_2(τ) * (3b + 2m * τ) - exp(-τ) * (3b + 2m * (τ + 1)))
end

"""
    exponential_integral_2(x)

Approximate second order exponential integral, E_2(x).  This stitches together several series
expansions to get an approximation which is accurate within 1% for all `x`
"""
function exponential_integral_2(x)  # this implementation could definitely be improved
    if x == 0
        1.0
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
    1 +
    ((log(x) + ℇ - 1) +
     (-0.5 + (0.08333333333333333 + (-0.013888888888888888 +
                                     0.0020833333333333333 * x) * x) * x) * x) * x
end
function _expint_large(x)
    invx = 1 / x
    exp(-x) * (1 + (-2 + (6 + (-24 + 120 * invx) * invx) * invx) * invx) * invx
end
function _expint_2(x)
    x -= 2
    0.037534261820486914 +
    (-0.04890051070806112 +
     (0.033833820809153176 +
      (-0.016916910404576574 +
       (0.007048712668573576 - 0.0026785108140579598 * x) * x) * x) * x) * x
end
function _expint_3(x)
    x -= 3
    0.010641925085272673 +
    (-0.013048381094197039 +
     (0.008297844727977323 +
      (-0.003687930990212144 + (0.0013061422257001345 - 0.0003995258572729822 * x) * x) * x) * x) *
    x
end
function _expint_4(x)
    x -= 4
    0.0031982292493385146 +
    (-0.0037793524098489054 +
     (0.0022894548610917728 +
      (-0.0009539395254549051 + (0.00031003034577284415 - 8.466213288412284e-5 * x) * x) * x) * x) *
    x
end
function _expint_5(x)
    x -= 5
    0.000996469042708825 +
    (-0.0011482955912753257 +
     (0.0006737946999085467 +
      (-0.00026951787996341863 + (8.310134632205409e-5 - 2.1202073223788938e-5 * x) * x) * x) * x) *
    x
end
function _expint_6(x)
    x -= 6
    0.0003182574636904001 +
    (-0.0003600824521626587 +
     (0.00020656268138886323 +
      (-8.032993165122457e-5 + (2.390771775334065e-5 - 5.8334831318151185e-6 * x) * x) * x) * x) * x
end
function _expint_7(x)
    x -= 7
    0.00010350984428214624 +
    (-0.00011548173161033826 +
     (6.513442611103688e-5 +
      (-2.4813114708966427e-5 + (7.200234178941151e-6 - 1.7027366981408086e-6 * x) * x) * x) * x) *
    x
end
function _expint_8(x)
    x -= 8
    3.413764515111217e-5 +
    (-3.76656228439249e-5 +
     (2.096641424390699e-5 +
      (-7.862405341465122e-6 + (2.2386015208338193e-6 - 5.173353514609864e-7 * x) * x) * x) * x) * x
end

end #module
