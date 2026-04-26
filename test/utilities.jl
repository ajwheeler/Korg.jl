using Printf: @sprintf

"""
    assert_allclose(actual, reference; [rtol], [atol], [err_msg], [error_location_fmt])

Raises an assertion exception when any entries in `actual` and `reference` are not equal up to the
desired tolerance. In other words, all entries of `actual` and `reference` must satisfy:
`abs(actual - reference) <= rtol * abs(reference) + atol`

This is inspired by numpy's assert_allclose function. However, we've changed some default values.
Additionally, in cases where actual and reference don't have equal values but do have the same
shape, the locations of the max errors are also printed.

# Keywords

  - `rtol`: The desired relative tolerance.  Use 0 to skip.
  - `atol`: The desired absolute tolerance.  Skipped by default.
  - `err_msg`: An optional error message that can be printed when the arrays are not close
  - `error_location_fmt`: An optional keyword used to specify a function that creates a string
    describing the location in the arrays where values are unequal, given a `CartesianIndex` object.
    This should not be specified in cases where actual and reference have different shapes.
  - `print_rachet_info`: set this to false to avoid printing a notice and stacktrace with rtol
    or atol can be tightened.
"""
function assert_allclose(actual, reference; rtol=1e-7, atol=0.0, err_msg=nothing,
                         error_location_fmt=nothing, print_rachet_info=true)
    # log if the comparison can be racheted down
    diff = abs.(actual .- reference)
    relative_diff = diff ./ abs.(reference)
    # handle the case where reference and computed value are both 0, avoid NaNs
    relative_diff[ForwardDiff.value.(diff).==0] .= 0

    if print_rachet_info
        if all(diff .< 0.1 * atol)
            @info "test can be racheted down: atol=$(atol), but the max diff is $(maximum(diff))"
            display(stacktrace())
        end
        if all(relative_diff .< 0.1 * rtol)
            @info ("test can be racheted down: rtol=$(rtol), but the max relative diff is "
                   *
                   "$(maximum(relative_diff))")
            display(stacktrace())
        end
    end

    # return promptly if the check passes
    if all((diff .== 0) .| (diff .<= rtol .* abs.(reference) .+ atol))
        return true
    end

    common_layout = size(actual) == size(reference)
    if isnothing(error_location_fmt)
        error_location_fmt = (arg) -> string((arg isa Integer) ? arg : arg.I)
    elseif !common_layout
        throw("error_location_fmt must be nothing b/c actual and reference have different shapes")
    end

    # now let's provide a detailed error message (this can take longer):

    # determine what the max error is (and possibly where it happended)
    diffmax = argmax(diff) #this in an index
    relmax = argmax(relative_diff) #this is an index

    # format the message
    lines = isnothing(err_msg) ? [] : [err_msg]
    push!(lines, "The arrays aren't consistent to within atol = $atol, rtol = $rtol")

    if (rtol != 0.0) || (rtol == atol == 0.0)
        err = @sprintf("Max Rel Diff:  %g = |%g - %g|/|%g|",
                       relative_diff[relmax], actual[relmax], reference[relmax], reference[relmax])
        err *= (common_layout) ? " at " * error_location_fmt(relmax) : ""
        push!(lines, err)
    end

    if (atol != 0.0) || (rtol == atol == 0.0)
        err = @sprintf("Max Abs Diff: %g = |%g - %g|",
                       diff[diffmax], actual[diffmax], reference[diffmax])
        err *= (common_layout) ? " at " * error_location_fmt(diffmax) : ""
        push!(lines, err)
    end

    throw(ErrorException(join(lines, "\n")))
end

"""
    assert_allclose_grid(actual, reference, independent_vars; kwargs...)

This function is a convenience wrapper around `assert_allclose` that provides slightly nicer error
messages based on the independent variables generating the grid.  It requires that `actual` and
`reference` have the same shape.

`independent_vars` should be a `ndims(actual)` element array. The ith element of `independent_vars`
should be 2-tuple or 3-tuple:

  - `independent_vars[i][1]` should specify the name of the independent variable along the ith axis
  - `independent_vars[i][2]` should specify the values of the independent variable along the ith axis
  - `independent_vars[i][3]`, if present, should specify the units of the independent variable.

For keyword arguments, see [`assert_allclose`](@ref).
"""
function assert_allclose_grid(actual, reference, independent_vars; kwargs...)
    @assert length(independent_vars) == ndims(actual)

    function error_location_fmt(index)
        ax_strs = map(1:ndims(actual)) do ax
            if length(independent_vars[ax]) ∉ [2, 3]
                error("independent_vars[$ax] should only have 2 or 3 elements")
            end

            name, vals = independent_vars[ax][1:2]
            unit_suffix = (length(independent_vars[ax]) == 3) ? " $(independent_vars[ax][3])" : ""

            if length(vals) != size(actual)[ax]
                error("length(independent_vars[$ax][2]) should be $(length(vals))")
            end

            "$name = $(vals[index[ax]])" * unit_suffix
        end
        join(ax_strs, ", ")
    end
    assert_allclose(actual, reference; error_location_fmt=error_location_fmt, kwargs...)
end

"""
Same as the other assert_allclose testing utils, but for dicts.  Raises an
assertion exception when any values in `actual` and `reference` are not equal up
to the desired tolerance or when their keys are not the same.

# Keywords

  - `rtol`: The desired relative tolerance. Use 0 to skip.
  - `atol`: The desired absolute tolerance. Skipped by default.
  - `err_msg`: An optional error message that can be printed when the dicts are not close.
  - `print_rachet_info`: set this to false to avoid printing a notice and stacktrace when rtol
    or atol can be tightened.
"""
function assert_allclose_dict(actual::AbstractDict, reference::AbstractDict;
                              rtol=1e-7, atol=0.0, err_msg=nothing, print_rachet_info=true)
    actual_keys = Set(keys(actual))
    reference_keys = Set(keys(reference))

    if actual_keys != reference_keys
        missing_keys = setdiff(reference_keys, actual_keys)
        extra_keys = setdiff(actual_keys, reference_keys)
        lines = isnothing(err_msg) ? String[] : [err_msg]
        push!(lines, "The dictionaries do not contain the same keys.")
        isempty(missing_keys) || push!(lines, "Missing keys: $(collect(missing_keys))")
        isempty(extra_keys) || push!(lines, "Extra keys: $(collect(extra_keys))")
        throw(ErrorException(join(lines, "\n")))
    end

    ks = collect(reference_keys)

    diff = Dict(k => abs(actual[k] - reference[k]) for k in ks)

    relative_diff = Dict{eltype(ks),Float64}()
    for k in ks
        if ForwardDiff.value(diff[k]) == 0
            relative_diff[k] = 0.0
        else
            relative_diff[k] = diff[k] / abs(reference[k])
        end
    end

    if print_rachet_info
        if all(values(diff) .< 0.1 * atol)
            @info "test can be racheted down: atol=$(atol), but the max diff is $(maximum(values(diff)))"
            display(stacktrace())
        end

        if all(values(relative_diff) .< 0.1 * rtol)
            @info("test can be racheted down: rtol=$(rtol), but the max relative diff is "*
                  "$(maximum(values(relative_diff)))")
            display(stacktrace())
        end
    end

    if all(k -> (diff[k] == 0) ||
               (diff[k] <= rtol * abs(reference[k]) + atol),
           ks)
        return true
    end

    diffmax = argmax(diff)
    relmax = argmax(relative_diff)

    lines = isnothing(err_msg) ? String[] : [err_msg]
    push!(lines, "The dictionaries aren't consistent to within atol = $atol, rtol = $rtol")

    if (rtol != 0.0) || (rtol == atol == 0.0)
        err = @sprintf("Max Rel Diff:  %g = |%g - %g|/|%g|",
                       relative_diff[relmax],
                       actual[relmax],
                       reference[relmax],
                       reference[relmax])
        err *= " at key $(repr(relmax))"
        push!(lines, err)
    end

    if (atol != 0.0) || (rtol == atol == 0.0)
        err = @sprintf("Max Abs Diff: %g = |%g - %g|",
                       diff[diffmax],
                       actual[diffmax],
                       reference[diffmax])
        err *= " at key $(repr(diffmax))"
        push!(lines, err)
    end

    throw(ErrorException(join(lines, "\n")))
end