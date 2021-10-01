using Korg


struct Interval # represents an exclusive interval
    lower::Union{Float64, Nothing} # nothing is equivalent to ∞
    upper::Union{Float64, Nothing} # nothing is equivalent to ∞

    function Interval(lower::Union{Float64, Nothing}, upper::Union{Float64, Nothing})
        if isinf(lower)
            @assert lower < 0 "The lower bound can't be +∞"
            lower = nothing
        end

        if isinf(upper)
            @assert upper > 0 "The upper bound can't be -∞"
            upper = nothing
        end

        if (!isnothing(lower)) && (!isnothing(upper)) && (lower >= upper)
            error("the upper bound must exceed the lower bound")
        end
        new(lower, upper)
    end
end

# TODO: re-add the ability to initialize from tuple

# take some care to avoid subnormal values (they can slow things down a lot
nextfloat_skipsubnorm(v::F) where {F<:AbstractFloat} =
    ifelse(-floatmin(F) ≤ v < 0, F(0), ifelse(0 ≤ v < floatmin(F), floatmin(F), nextfloat(v)))
prevfloat_skipsubnorm(v::F) where {F<:AbstractFloat} =
    ifelse(-floatmin(F) < v ≤ 0, -floatmin(F), ifelse(0 < v ≤ floatmin(F), F(0), prevfloat(v)))

clamp(v, b::Interval) = ifelse(v<b.lower, nextfloat_skipsubnorm(v),
                               ifelse(v > b.upper, prevfloat_skipsubnorm(v), v))

# we should probably drop this
function Interval(str::AbstractString)
    # the regex expression is defined in "comment mode" (the expression ignores all whitespace
    # not escaped with `\` and any text between a `#` and a newline)
    re = r"(?<lbracket>[\[(])\s*                                              # opening bracket
           (?<lower>([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)|([-+]?∞)),\s*   # left number
           (?<upper>([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)|([-+]?∞))\s*    # right number
           (?<rbracket>[\])])"x                                               # closing bracket
    m = match(re, str)

    inf_dict = Dict("-∞"=> -Inf, "∞"=> Inf, "+∞"=> Inf)

    lower = (m[:lower] in keys(inf_dict)) ? inf_dict[m[:lower]] : parse(Float64, m[:lower])
    if isfinite(lower) && m[:lbracket] == "[" # otherwise, m[:lbracket] is equal to "("
        lower = prevfloat_skipsubnorm(lower)
    end

    upper = (m[:upper] in keys(inf_dict)) ? inf_dict[m[:upper]] : parse(Float64, m[:upper])
    if isfinite(lower) && m[:rbracket] == "]" # otherwise, m[:rbracket] is equal to ")"
        upper = nextfloat_skipsubnorm(upper)
    end

    Interval(lower, upper)
end


function contained(value::Real, interval::Interval)
    ( (isnothing(interval.lower) || (value > interval.lower)) &&
      (isnothing(interval.upper) || (value < interval.upper)) &&
      isfinite(value) ) # this can give the wrong result if isfinite(value) is omitted
end

function contained_slice(monotonic_vals::AbstractRange, interval::Interval)
    error("this needs to be tested")
    if monotonic_vals[1] < monotonic_vals[2]
        # monotonically increasing
        start = searchsortedfirst(interval.lower)
        stop = searchsortedlast(interval.upper)
    else
        start = searchsortedlast(interval.upper; rev = true)
        stop = searchsortedfirst(interval.lower; rev = true)
    end
    start, stop
end

function _convert_λ_to_ν_bound(λ_bound::Interval)
    calc_λ(ν) = c_cgs/ν

    ν_upper = (λ_bound.lower == 0) ? nothing : Korg.c_cgs/λ_bound.lower
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

    ν_lower = isnothing(λ_bound.upper) ? 0.0 : Korg.c_cgs/λ_bound.upper
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

"""
    bounds_checked_absorption(func, photon_bound, photon_bound_type, temperature_bound)

Constructs a wrapped function that implements bounds checking and extrapolation

# Parameters
- `func`: a function that has a signature `f(ν::Real, T::Real, args...)::Real`, where `ν` is 
  frequency (in Hz) and `T` is temperature (in K)
- `photon_bound::Interval`: Interval of photon properties over which `func` is valid
- `λ_based_photon_bound::AbstractString`: This can be "λ", "ν", or "ν/T" which respectively mean
  that photon_bound corresponds to wavelength (in cm), frequency (Hz), or frequency divided by
  temperature (in Hz per K)
- `temperature_bound::Interval`: Interval of temperatures (in K) over which `func` is valid.

"""
function bounds_checked_absorption(func, photon_bound::Interval, photon_bound_type::AbstractString,
                                   temperature_bound::Interval)

    # if generating the inner function in this way involves too much overhead, we could:
    # - convert this outer function to a macro.
    # - perform the boundary checks inside of total_continuum_absorption

    ν_bound = if photon_bound_type == "λ"
        # It would probably be better (in terms of round-off errors) to convert ν to λ (instead of
        # converting the bounds) since that's the direction of the conversion within the function
        _convert_λ_to_ν_bound(photon_bound)
    elseif photon_bound_type == "ν"
        photon_bound
    elseif photon_bound_type != "ν/T"
        error(string("unrecognized photon_bound_type: ", photon_bound_type))
    else
        nothing
    end

    # If we want to support T to be a Vector:
    # - we'll have to assume that any vectors in `args...` are the same length as T.
    # - if we made ionization_energy a keyword argument, then you could assume that all entries in
    #   `args...` are vectors with the same length as T.
    # It's probably unnecessary to handle this case anyway, since the outputs should be contiguous
    # in ν to facillitate fast interpolation
    function wrapped_func(ν::AbstractVector{<:Real}, T::Real, args...;
                          extrapolate_bc = 0.0, out_α::Union{Nothing,AbstractVector} = nothing)

        if photon_bound_type == "ν/T"
            ν_bound = Interval(photon_bound.lower * T, photon_bound.upper * T)
        end

        # the following boundary condition handling is slow. It's mostly for debugging, anyways
        if isnothing(extrapolate_bc)
            bad_ν_index = findfirst(.!contained.(ν, Ref(ν_bound)))
            if !contained(T, temperature_bound)
                throw(DomainError(T, ("$func: invalid temperature. It should lie between "*
                                      "$(temperature_bound.lower) & $(temperature_bound.upper)")))
            elseif photon_bound_type == "ν/T" && !isnothing(bad_ν_index)
                throw(DomainError((ν[bad_ν_index], T),
                                  ("$func: invalid (ν, T) pair. ν/T should lie between "*
                                   "$(photon_bound.lower) & $(photon_bound.upper)")))
            elseif !isnothing(bad_ν_index)
                throw(DomainError(ν[bad_ν_index], ("$func: invalid freq. It should lie between "*
                                                   "$(ν_bound.lower) & $(ν_bound.upper)")))
            end
            extrapolate_bc = 0.0
        end

        if isnothing(out_α) # zero-initialize since we add computed values to out_α in place
            out_α = zeros(promote_type(eltype(ν), typeof(T)), length(ν))
        end

        if contained(T, temperature_bound)
            if ν isa AbstractRange
                first, last = contained_slice(ν, ν_bound)
                view(out_α, 1:(first-1)) .+= extrapolate_bc
                view(out_α, first:last) .+= func.(ν, T, args...)
                view(out_α, (last+1):length(out_α)) .+= extrapolate_bc
            else
                # if very specific cases (but probably not the general case), the commented code
                # might be faster (it depends on which function we're calling). clamp is essential,
                # because ifelse eagerly evaluates the outputs (unlike the ternary operator)
                #out_α .+= ifelse.(contained.(ν, Ref(ν_bound)),
                #                  func.(clamp.(ν, Ref(ν_bound)), T, args...),
                #                  extrapolate_bc)
                for i in 1:length(ν)
                    out_α[i] += contained(ν[i], ν_bound) ? func(ν[i], T, args...) : extrapolate_bc
                end
            end
        else
            out_α .+= extrapolate_bc
        end
        out_α
    end

    # handle the case where ν is a scalar to avoid breaking tests. My intention is to remove this
    # if we ultimately decide to go with this function wrapping approach
    function wrapped_func(ν::Real, T::Real, args...; extrapolate_bc = 0.0)
        tmp = zeros(promote_type(typeof(ν), typeof(T)), 1)
        wrapped_func([ν], T, args...; extrapolate_bc = extrapolate_bc)[1]
    end

    wrapped_func
end
