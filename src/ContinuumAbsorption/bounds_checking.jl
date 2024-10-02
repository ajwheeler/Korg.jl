function _prep_bound_err(ν, T, ν_indices, func, ν_bound, temperature_bound)
    signs = "-+"
    fmt(x) = isinf(x) ? "$(signs[1+(x>0)])∞" : "$x"

    if !contained(T, temperature_bound)
        DomainError(T,
                    ("$func: invalid temperature. It should lie between " *
                     "$(temperature_bound.lower) & $(temperature_bound.upper)"))
    else
        bad_ν = (1 in ν_indices) ? ν[last(ν_indices)+1] : ν[1]
        DomainError(bad_ν,
                    ("$func: invalid frequency. It should lie between " *
                     "$(fmt(ν_bound.lower)) Hz & $(fmt(ν_bound.upper)) Hz"))
    end
end

"""
    bounds_checked_absorption(func; ν_bound, temp_bound)

Constructs a wrapped function that implements bounds checking and extrapolation

# Parameters

  - `func`: a function that has a signature `f(ν::Real, T::Real, args...)::Real`, where `ν` is
    frequency (in Hz) and `T` is temperature (in K)
  - `ν_bound::Interval`: Interval of frequencies (in Hz) over which `func` is valid.
  - `temp_bound::Interval`: Interval of temperatures (in K) over which `func` is valid.

The resulting function will have this signature:

```
wrapped_func(ν::AbstractVector{<:Real}, T::Real, args...; kwargs...)
```

# Wrapped Function Parameters

  - `ν::AbstractVector{<:Real}`: sorted vector of frequencies (in Hz)
  - `T::Real`: temperature (in K)
  - `args...`: function-specific arguments

For a description for the `kwargs...`, see [Continuum Absorption Kwargs](@ref).    # if generating the inner function in this way involves too much overhead, we could:
"""
function bounds_checked_absorption(func; ν_bound=Interval(0, Inf), temp_bound=Interval(0, Inf))
    # if generating the inner function in this way involves too much overhead, we could:
    # - convert this outer function to a macro.
    # - perform the boundary checks inside of total_continuum_absorption

    # It's probably unnecessary to let T be a Vector since the outputs should be contiguous in ν to
    # fast interpolation. Broadcasting over T would be fine despite the fact that it entails
    # redundant bounds checking because real problems have only order ~100 T values
    function wrapped_func(ν::AbstractVector{<:Real}, T::Real, args...; error_oobounds::Bool=false,
                          out_α::Union{Nothing,AbstractVector}=nothing)
        @assert issorted(ν, rev=first(ν) > last(ν)) "ν should be sorted"

        α_type = promote_type(eltype(ν), typeof(T), typeof.(args)...)
        out_α = isnothing(out_α) ? zeros(α_type, length(ν)) : out_α
        @assert (eltype(out_α) == α_type) && (length(out_α) == length(ν))

        # find indices where out_α can be updated from arguments that are in bounds
        idx = (!contained(T, temp_bound)) ? (1:0) : contained_slice(ν, ν_bound)

        # handle any out-of-bounds arguments. If !error_oobounds, assume that there is 0 absorption
        # for these args. Otherwise, throw error
        if (idx != 1:length(ν)) && error_oobounds  # (it's okay for this to be slow)
            throw(_prep_bound_err(ν, T, idx, func, ν_bound, temp_bound))
        end

        view(out_α, idx) .+= func.(ν[idx], T, args...)
        out_α
    end

    wrapped_func
end
