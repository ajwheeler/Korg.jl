# take some care to avoid subnormal values (they can slow things down a lot!)
nextfloat_skipsubnorm(v::F) where {F<:AbstractFloat} =
    ifelse(-floatmin(F) ≤ v < 0, F(0), ifelse(0 ≤ v < floatmin(F), floatmin(F), nextfloat(v)))
prevfloat_skipsubnorm(v::F) where {F<:AbstractFloat} =
    ifelse(-floatmin(F) < v ≤ 0, -floatmin(F), ifelse(0 < v ≤ floatmin(F), F(0), prevfloat(v)))

struct Interval # represents an exclusive interval
    lower::Float64
    upper::Float64

    function Interval(lower::Real, upper::Real; exclusive_lower = true, exclusive_upper = true)
        @assert lower < upper "the upper bound must exceed the lower bound"
        lower, upper = Float64(lower), Float64(upper)
        new((exclusive_lower || isinf(lower)) ? lower : prevfloat_skipsubnorm(lower),
            (exclusive_upper || isinf(upper)) ? upper : nextfloat_skipsubnorm(upper))
    end
end

# convenience function for defining interval where both bounds are inclusive
closed_interval(lo, up) = Interval(lo, up; exclusive_lower = false, exclusive_upper = false)

"""
    contained(value, interval)

Returns whether `value` is contained by `interval`.
"""
function contained(value::Real, interval::Interval)
    # this can give the wrong result if isfinite(value) is omitted
    (value > interval.lower) && (value < interval.upper) && isfinite(value)
end


"""
    contained_slice(vals, interval)

Returns a range of indices denoting the elements of `vals` that are contained by `interval`. When
no entries are contained by interval, this returns `(1,0)` (which is a valid empty slice).
This function assumes that `vals` is sorted in increasing/decreasing order.
"""
function contained_slice(vals::AbstractVector, interval::Interval)
    if first(vals) < last(vals)
        first_ind = searchsortedfirst(vals, interval.lower)
        last_ind = searchsortedlast(vals, interval.upper)
    else
        first_ind = searchsortedlast(vals, interval.upper; rev = true) + 1
        last_ind = searchsortedfirst(vals, interval.lower; rev = true) - 1
    end
    first_ind:last_ind
end

# inbound_ν_neighbor
# oobound_ν_neighbor

function _λ_to_ν_bound(λ_bound::Interval)
    calc_λ(ν) = c_cgs/ν

    ν_upper = (λ_bound.lower == 0) ? Inf : c_cgs/λ_bound.lower
    if isfinite(ν_upper)
        # adjust ν_upper so that we can pass the assert statement
        while calc_λ(prevfloat(ν_upper)) ≤ λ_bound.lower # ν_upper is too large, make it smaller
            ν_upper = prevfloat(ν_upper)
        end
        while calc_λ(ν_upper) > λ_bound.lower # ν_upper is too small, make it larger
            ν_upper = nextfloat(ν_upper)
        end
        @assert (calc_λ(ν_upper) ≤ λ_bound.lower) && (calc_λ(prevfloat(ν_upper)) > λ_bound.lower)
    end

    ν_lower = isnothing(λ_bound.upper) ? 0.0 : c_cgs/λ_bound.upper
    if isfinite(ν_lower) && (!isnothing(λ_bound.upper)) && isfinite(λ_bound.upper)
        # adjust ν_lower so that we can pass the assert statement
        while calc_λ(nextfloat(ν_lower)) ≥ λ_bound.upper # ν_lower is too small, make it larger
            ν_lower = nextfloat(ν_lower)
        end
        while calc_λ(ν_lower) < λ_bound.upper # ν_lower is too large, make it smaller
            ν_lower = prevfloat(ν_lower)
        end
        @assert (calc_λ(ν_lower) ≥ λ_bound.upper) && (calc_λ(nextfloat(ν_lower)) < λ_bound.upper)
    end

    Interval(ν_lower, ν_upper)
end

function _prep_bound_err(ν, T, ν_indices, func, ν_bound, ν_div_T_bound, temperature_bound)
    signs = "-+"
    fmt(x) = isinf(x) ? "$(signs[1+(x>0)])∞" : "$x"

    if !contained(T, temperature_bound)
        DomainError(T, ("$func: invalid temperature. It should lie between "*
                        "$(temperature_bound.lower) & $(temperature_bound.upper)"))
    else
        bad_ν = (1 in ν_indices) ? ν[last(ν_indices) + 1] : ν[1]
        if !isnothing(ν_div_T_bound)
            DomainError((bad_ν, T), ("$func: invalid (ν, T) pair. ν/T should lie between "*
                                     "$(ν_div_T_bound.lower) & $(ν_div_T_bound.upper)"))
        else
            DomainError(bad_ν, ("$func: invalid freq. It should lie between "*
                                "$(fmt(ν_bound.lower)) Hz & $(fmt(ν_bound.upper)) Hz"))
        end
    end
end


"""
    bounds_checked_absorption(func; ν_bound, ν_div_T_bound, temp_bound)

Constructs a wrapped function that implements bounds checking and extrapolation

# Parameters
- `func`: a function that has a signature `f(ν::Real, T::Real, args...)::Real`, where `ν` is 
  frequency (in Hz) and `T` is temperature (in K)
- `photon_bound::Interval`: Interval of photon properties over which `func` is valid
- `λ_based_photon_bound::AbstractString`: This can be "λ", "ν", or "ν/T" which respectively mean
  that photon_bound corresponds to wavelength (in cm), frequency (Hz), or frequency divided by
  temperature (in Hz per K)
- `temp_bound::Interval`: Interval of temperatures (in K) over which `func` is valid.


The resulting function will have a signature that loosely reflects
wrapped_func(ν::AbstractVector{<:Real}, T::Real, args...;
             extrapolate_bc = 0.0, out_α::Union{Nothing, AbstractVector{<:Real}})

# Wrapped Parameters
- `ν::AbstractVector{<:Real}`: sorted vector (either forward or reverse order) of frequencies

"""
function bounds_checked_absorption(func; ν_bound::Union{Interval,Nothing} = nothing,
                                   ν_div_T_bound::Union{Interval,Nothing} = nothing,
                                   temp_bound::Union{Interval,Nothing} = nothing)

    # if generating the inner function in this way involves too much overhead, we could:
    # - convert this outer function to a macro.
    # - perform the boundary checks inside of total_continuum_absorption

    if isnothing(ν_bound) && isnothing(ν_div_T_bound)
        error("Either ν_bound or ν_div_T_bound must be specified")
    elseif (!isnothing(ν_bound)) && (!isnothing(ν_div_T_bound))
        error("Current implementation doesn't allow both ν_bound & ν_div_T_bound to be specified")
    end
    @assert !isnothing(temp_bound) "the current implementation requires temp_bound"

    # It's unnecessary to let T be a Vector since the outputs should be contiguous in ν to
    # fast interpolation. But, if we ever wanted to supporth that case:
    # - we'll have to assume that any vectors in `args...` are the same length as T.
    # - if we made ionization_energy a keyword argument, then you could assume that all entries in
    #   `args...` are vectors with the same length as T.
    function wrapped_func(ν::AbstractVector{<:Real}, T::Real, args...;
                          extrapolate_bc = 0.0, out_α::Union{Nothing,AbstractVector} = nothing)
        @assert issorted(ν, rev = first(ν) > last(ν)) "ν should be sorted"

        α_type = reduce(promote_type, [eltype(ν), typeof(T), mapreduce(typeof, promote_type, args)])
        out_α = isnothing(out_α) ? zeros(α_type, length(ν)) : out_α
        @assert (eltype(out_α) == α_type) && (length(out_α) == length(ν))

        # find indices where out_α can be updated from arguments that are in bounds
        tmp_ν_bound = isnothing(ν_div_T_bound) ? ν_bound :
            Interval(ν_div_T_bound.lower * T, ν_div_T_bound.upper * T)
        idx = (!contained(T, temp_bound)) ? (1:0) : contained_slice(ν, tmp_ν_bound)

        if length(idx) != length(ν) # handle arguments that are out of bounds
            if isnothing(extrapolate_bc) # throw error, (it's okay for this to be slow)
                throw(_prep_bound_err(ν, T, idx, func, ν_bound, ν_div_T_bound, temp_bound))
            else
                view(out_α, 1:(first(idx)-1)) .+= extrapolate_bc
                view(out_α, (last(idx)+1):length(ν)) .+= extrapolate_bc
            end
        end

        view(out_α, idx) .+= func.(ν[idx], T, args...)
        out_α
    end

    wrapped_func
end
