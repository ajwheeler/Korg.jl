using Statistics: quantile
using Interpolations: linear_interpolation, Flat
using SparseArrays: spzeros

normal_pdf(Δ, σ) = exp(-0.5 * Δ^2 / σ^2) / √(2π) / σ

"""
    merge_bounds(bounds, merge_distance=0.0)

Sort a vector of lower-bound, upper-bound pairs and merge overlapping ranges.

Returns a pair containing:

  - a vector of merged bounds
  - a vector of vectors of indices of the original bounds which were merged into each merged bound
"""
function merge_bounds(bounds, merge_distance=0.0)
    bound_indices = 1:length(bounds)

    # sort by lower bound
    s = sortperm(bounds; by=first)
    bounds = bounds[s]
    bound_indices = bound_indices[s]

    new_bounds = [bounds[1]]
    indices = [[bound_indices[1]]]
    for i in 2:length(bounds)
        # if these bounds are within merge_distance of the previous, extend the previous,
        # otherwise add them to the list
        if bounds[i][1] <= new_bounds[end][2] + merge_distance
            new_bounds[end] = (new_bounds[end][1], max(bounds[i][2], new_bounds[end][2]))
            push!(indices[end], bound_indices[i])
        else
            push!(new_bounds, bounds[i])
            push!(indices, [bound_indices[i]])
        end
    end
    new_bounds, indices
end

# Convert R to a value based on its type
# used in `_lsf_bounds_and_kernel`
_resolve_R(R::Real, λ0) = R
_resolve_R(R::Function, λ0) = R(λ0 * 1e8)  # R is a function of λ in Å

# Core LSF calculation shared by all variants
function _lsf_bounds_and_kernel(synth_wls::Wavelengths, λ0, R, window_size)
    R_val = _resolve_R(R, λ0)
    σ = λ0 / R_val / (2sqrt(2log(2))) # convert Δλ = λ0/R (FWHM) to sigma

    # Calculate bounds and kernel
    lb = searchsortedfirst(synth_wls, λ0 - window_size * σ)
    ub = searchsortedlast(synth_wls, λ0 + window_size * σ)
    @views ϕ = normal_pdf.(synth_wls[lb:ub] .- λ0, σ)
    normalized_ϕ = ϕ ./ sum(ϕ)

    lb, ub, normalized_ϕ
end

"""
    apply_LSF(flux, wls, R; window_size=4)

Applies a gaussian line spread function the the spectrum with flux vector `flux` and wavelengths
`wls` with constant spectral
resolution (, ``R = \\lambda/\\Delta\\lambda``, where ``\\Delta\\lambda`` is the LSF FWHM.  The
`window_size` argument specifies how far out to extend the convolution kernel in standard deviations.

For the best match to data, your wavelength range should extend a couple ``\\Delta\\lambda`` outside
the region you are going to compare.

If you are convolving many spectra defined on the same wavelenths to observational resolution, you
will get much better performance using [`compute_LSF_matrix`](@ref).

# Arguments

  - `flux`: the flux vector to convolve
  - `wls`: wavelengths in any format [described here](@ref wldocs)
  - `R`: the resolving power, ``R = \\lambda/\\Delta\\lambda``

# Keyword Arguments

  - `window_size` (default: 4): how far out to extend the convolution kernel in units of the LSF width (σ, not HWHM)
  - `renormalize_edge` (default: `true`): whether or not to renormalize the LSF at the edge of the wl
    range.  This doen't matter as long as `synth_wls` extends to large and small enough wavelengths.

!!! warning

      - This is a naive, slow implementation.  Do not use it when performance matters.

      - `apply_LSF` will have weird behavior if your wavelength grid is not locally linearly-spaced.
        It is intended to be run on a fine wavelength grid, then downsampled to the observational (or
        otherwise desired) grid.
"""
function apply_LSF(flux::AbstractVector{F}, wls, R; window_size=4) where F<:Real
    if R == Inf
        return copy(flux)
    end
    wls = Wavelengths(wls)
    convF = zeros(F, length(flux))
    for i in eachindex(wls)
        λ0 = wls[i]
        lb, ub, normalized_ϕ = _lsf_bounds_and_kernel(wls, λ0, R, window_size)
        convF[i] = sum(flux[lb:ub] .* normalized_ϕ)
    end
    convF
end
@deprecate constant_R_LSF(args...; kwargs...) apply_LSF(args...; kwargs...)

"""
    compute_LSF_matrix(synth_wls, obs_wls, R; kwargs...)

Construct a sparse matrix, which when multiplied with a flux vector defined over wavelenths
`synth_wls`, applies a gaussian line spead function (LSF) and resamples to the wavelenths `obs_wls`.

# Arguments

  - `synth_wls`: the synthesis wavelengths in any form [described here](@ref wldocs)
  - `obs_wls`: the wavelengths of the observed spectrum
  - `R`: the resolving power, ``R = \\lambda/\\Delta\\lambda``

# Keyword Arguments

  - `window_size` (default: 4): how far out to extend the convolution kernel in units of the LSF width (σ, not HWHM)
  - `verbose` (default: `true`): whether or not to emit warnings and information to stdout/stderr.
  - `step_tolerance`: the maximum difference between adjacent wavelengths in `synth_wls` for them to be
    considered linearly spaced.  This is only used if `synth_wls` is a vector of wavelengths rather
    than a range or vector or ranges.

For the best match to data, your wavelength range should extend a couple ``\\Delta\\lambda`` outside
the region you are going to compare.

[`Korg.apply_LSF`](@ref) can apply an LSF to a single flux vector efficiently. This function is
relatively slow, but once the LSF matrix is constructed, convolving spectra to observational
resolution via matrix multiplication is fast.
"""
function compute_LSF_matrix(synth_wls, obs_wls, R; window_size=4, verbose=true)
    if first(obs_wls) >= 1
        # don't change obs_wls in place, make a copy
        obs_wls = obs_wls / 1e8 # Å to cm
    end
    synth_wls = Wavelengths(synth_wls)
    if verbose &&
       !((first(synth_wls) - 0.01) <= first(obs_wls) <= last(obs_wls) <= (last(synth_wls) + 0.01))
        @warn "Synthesis wavelenths $(synth_wls) are not superset of observation wavelenths" *
              " ($(first(obs_wls)*1e8) Å—$(last(obs_wls)*1e8) Å) in LSF matrix."
    end
    LSF = spzeros((length(synth_wls), length(obs_wls)))
    for i in eachindex(obs_wls)
        λ0 = obs_wls[i]
        lb, ub, normalized_ϕ = _lsf_bounds_and_kernel(synth_wls, λ0, R, window_size)
        LSF[lb:ub, i] .+= normalized_ϕ
    end
    LSF'
end

"""
    apply_rotation(flux, wls, vsini, ε=0.6)

Given a spectrum `flux` sampled at wavelengths `wls` for a non-rotating star, compute the spectrum
that would emerge given projected rotational velocity `vsini` (in km/s) and linear limb-darkening
coefficient `ε`: ``I(\\mu) = I(1) (1 - \\varepsilon + \varepsilon \\mu))``.  See, for example,
Gray equation 18.14.

# Arguments

  - `flux`: the flux vector to rotate
  - `wls`: wavelengths in any format [described here](@ref wldocs)
  - `vsini`: projected rotational velocity in km/s
  - `ε`: linear limb-darkening coefficient (default: 0.6)
"""
function apply_rotation(flux, wls, vsini, ε=0.6)
    wls = Wavelengths(wls)
    newflux = similar(flux)
    lower_index = 1
    upper_index = length(wls.wl_ranges[1])
    newflux[lower_index:upper_index] .= _apply_rotation_core(view(flux, lower_index:upper_index),
                                                             wls.wl_ranges[1], vsini, ε)
    for i in 2:length(wls.wl_ranges)
        lower_index = upper_index + 1
        upper_index = lower_index + length(wls.wl_ranges[i]) - 1
        newflux[lower_index:upper_index] .= _apply_rotation_core(view(flux,
                                                                      lower_index:upper_index),
                                                                 wls.wl_ranges[i], vsini, ε)
    end
    newflux
end
function _apply_rotation_core(flux, wls::StepRangeLen, vsini, ε=0.6)
    if vsini == 0
        return copy(flux)
    end

    if first(wls) > 1
        wls *= 1e-8 # Å to cm
    end
    vsini *= 1e5 # km/s to cm/s

    newFtype = promote_type(eltype(flux), eltype(wls), typeof(vsini), typeof(ε))
    newF = zeros(newFtype, length(flux))

    # precompute constants
    c1 = 2(1 - ε)
    c2 = π * ε / 2
    # c3 is the denomicator.  the factor of v_L in Gray becomes Δλrot (because we are working in
    # wavelenths) and moves inside the loop
    c3 = π * (1 - ε / 3)

    for i in 1:length(flux)
        Δλrot = wls[i] * vsini / Korg.c_cgs # Å

        lb = searchsortedfirst(wls, wls[i] - Δλrot)
        ub = searchsortedlast(wls, wls[i] + Δλrot)
        Fwindow = flux[lb:ub]

        detunings = [-Δλrot; (lb-i+1/2:ub-i-1/2) * step(wls); Δλrot]

        ks = _rotation_kernel_integral_kernel.(c1, c2, c3 * Δλrot, detunings, Δλrot)
        newF[i] = sum(ks[2:end] .* Fwindow) - sum(ks[1:end-1] .* Fwindow)
    end
    newF
end

# the indefinite integral of the rotation kernel
function _rotation_kernel_integral_kernel(c1, c2, c3, detuning, Δλrot)
    if abs(detuning) == Δλrot
        return sign(detuning) * 0.5 # make it nan-safe for ForwardDiff
    end
    (0.5 * c1 * detuning * sqrt(1 - detuning^2 / Δλrot^2)
     + 0.5 * c1 * Δλrot * asin(detuning / Δλrot)
     + c2 * (detuning - detuning^3 / (3 * Δλrot^2))) / c3
end

"""
    air_to_vacuum(λ; cgs=λ<1)

Convert λ from an air to vacuum.  λ is assumed to be in Å if it is ⩾ 1, in cm otherwise.  Formula
from Birch and Downs (1994) via the VALD website.

See also: [`vacuum_to_air`](@ref).
"""
function air_to_vacuum(λ; cgs=λ < 1)
    if cgs
        λ *= 1e8
    end
    s = 1e4 / λ
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s^2) +
        0.0001599740894897 / (38.92568793293 - s^2)
    if cgs
        λ *= 1e-8
    end
    λ * n
end

"""
    vacuum_to_air(λ; cgs=λ<1)

convert λ from a vacuum to air.  λ is assumed to be in Å if it is ⩾ 1, in cm otherwise.  Formula
from Birch and Downs (1994) via the VALD website.

See also: [`air_to_vacuum`](@ref).
"""
function vacuum_to_air(λ; cgs=λ < 1)
    if cgs
        λ *= 1e8
    end
    s = 1e4 / λ
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s^2) + 0.00015998 / (38.9 - s^2)
    if cgs
        λ *= 1e-8
    end
    λ / n
end

# take some care to avoid subnormal values (they can slow things down a lot!)
function _nextfloat_skipsubnorm(v::F) where {F<:AbstractFloat}
    ifelse(-floatmin(F) ≤ v < 0, F(0), ifelse(0 ≤ v < floatmin(F), floatmin(F), nextfloat(v)))
end
function _prevfloat_skipsubnorm(v::F) where {F<:AbstractFloat}
    ifelse(-floatmin(F) < v ≤ 0, -floatmin(F), ifelse(0 < v ≤ floatmin(F), F(0), prevfloat(v)))
end

struct Interval # represents an exclusive interval
    lower::Float64
    upper::Float64

    function Interval(lower::Real, upper::Real; exclusive_lower=true, exclusive_upper=true)
        @assert lower<upper "the upper bound must exceed the lower bound"
        lower, upper = Float64(lower), Float64(upper)
        new((exclusive_lower || isinf(lower)) ? lower : _prevfloat_skipsubnorm(lower),
            (exclusive_upper || isinf(upper)) ? upper : _nextfloat_skipsubnorm(upper))
    end
end

# convenience function for defining interval where both bounds are inclusive
closed_interval(lo, up) = Interval(lo, up; exclusive_lower=false, exclusive_upper=false)

"""
    contained(value, interval)

Returns whether `value` is contained by `interval`.

# Examples

```julia-repl
julia> contained(0.5, Interval(1.0, 10.0))
false

julia> contained(5.0, Interval(1.0, 10.0))
true
```
"""
contained(value::Real, interval::Interval) = interval.lower < value < interval.upper

"""
    contained_slice(vals, interval)

Returns a range of indices denoting the elements of `vals` (which are assumed to be sorted in
increasing order) that are contained by `interval`. When no entries are contained by interval,
this returns `(1,0)` (which is a valid empty slice).
"""
contained_slice(vals::AbstractVector, interval::Interval) = searchsortedfirst(vals, interval.lower):searchsortedlast(vals,
                                                                                                                     interval.upper)

function _convert_λ_endpoint(λ_endpoint::AbstractFloat, λ_lower_bound::Bool)
    # determine the functions that:
    # - retrieve the neighboring ν val in the in-bounds and out-of-bounds directions
    # - specify the desired relationship between an in-bounds λ and λ_lower_bound
    inbound_ν_neighbor, oobound_ν_neighbor, inbound_λ_to_endpoint_relation = if λ_lower_bound
        (prevfloat, nextfloat, >)
    else
        (nextfloat, prevfloat, <)
    end

    ν_endpoint = (λ_endpoint == 0) ? Inf : c_cgs / λ_endpoint
    if isfinite(ν_endpoint) && (ν_endpoint != 0)
        # Adjust ν_endpoint in 2 cases:
        # Case 1: The neighboring ν to ν_endpoint that should be in-bounds, corresponds to a λ
        # value that is out-of-bounds. Nudge ν_endpoint in the "in-bounds" direction until resolved
        while !inbound_λ_to_endpoint_relation(c_cgs / inbound_ν_neighbor(ν_endpoint), λ_endpoint)
            ν_endpoint = inbound_ν_neighbor(ν_endpoint)
        end
        # Case 2: The current value of ν_endpoint corresponds to a λ that is in-bounds. Nudge
        # ν_endpoint in the "out-of-bounds" direction until resolved.
        while inbound_λ_to_endpoint_relation(c_cgs / ν_endpoint, λ_endpoint)
            ν_endpoint = oobound_ν_neighbor(ν_endpoint)
        end
    end
    ν_endpoint
end

"""
    λ_to_ν_bound(λ_bound)

Converts a λ `Inverval` (in cm) to an equivalent ν `Interval` (in Hz), correctly accounting for
tricky floating point details at the bounds.
"""
λ_to_ν_bound(λ_bound::Interval) = Interval(_convert_λ_endpoint(λ_bound.upper, false),
                                           _convert_λ_endpoint(λ_bound.lower, true))
