"""
    assert_allclose(actual, reference; [rtol], [atol], [err_msg], [error_location_fmt])

Raises an assertion exception when any entries in `actual` and `reference` are not equal up to the
desired tolerance. In other words, all entries of `actual` and `reference` must satisfy:
`abs(actual - reference) <= rtol * abs(reference) + atol`

This is inspired by numpy's assert_allclose function. However, we've changed some default values.
Additionally, in cases where actual and reference don't have equal values but do have the same 
shape, the locations of the max errors are also printed.

# Keywords
- `rtol`: The desired relative tolerance.  Use 0 to skip relative comparison.
- `atol`: The desired absolute tolerance
- `err_msg`: An optional error message that can be printed when the arrays are not close
- `error_location_fmt`: An optional keyword used to specify a function that creates a string 
  describing the location in the arrays where values are unequal, given a `CartesianIndex` object.
  This should not be specified in cases where actual and reference have different shapes.
"""
function assert_allclose(actual, reference; rtol = 1e-7, atol = 0.0, err_msg = nothing,
                         error_location_fmt = nothing)
    # make the initial check fast!
    if all( abs.(actual .- reference) .<= (rtol .* abs.(reference) .+ atol))
        return true
    end

    common_layout = size(actual) == size(reference)
    if isnothing(error_location_fmt)
        error_location_fmt = (arg) -> string( (arg isa Integer) ? arg : arg.I)
    elseif !common_layout
        throw("error_location_fmt must be nothing b/c actual and reference have different shapes")
    end

    # now let's provide a detailed error message (this can take longer):

    # determine what the max error is (and possibly where it happended)
    diff = abs.(actual .- reference)
    max_abs_diff_loc = argmax(diff)
    max_abs_diff = diff[max_abs_diff_loc]

    relative_diff = diff ./ abs.(reference)
    max_rel_error_magnitude_loc = argmax(relative_diff)
    max_rel_error_magnitude = relative_diff[max_rel_error_magnitude_loc]

    # format the message
    lines = isnothing(err_msg) ? [] : [err_msg]
    push!(lines, "The arrays aren't consistent to within atol = $atol, rtol = $rtol")

    if (rtol != 0.0) || (rtol == atol == 0.0)
        suffix = (common_layout) ? " at " * error_location_fmt(max_rel_error_magnitude_loc) : ""
        push!(lines, string("Max Rel Diff Magnitude: ",
                            isinf(max_rel_error_magnitude) ? "∞" : max_rel_error_magnitude,
                            suffix))
    end

    if (atol != 0.0) || (rtol == atol == 0.0)
        suffix = (common_layout) ? " at " * error_location_fmt(max_abs_diff_loc) : ""
        push!(lines, string("Max Abs Diff Magnitude: ", max_abs_diff, suffix))
    end

    throw(ErrorException(join(lines, "\n")))
    false
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
            if length(independent_vars[ax]) ∉ [2,3]
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
    assert_allclose(actual, reference; error_location_fmt = error_location_fmt, kwargs...)
end
