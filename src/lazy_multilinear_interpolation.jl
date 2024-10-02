struct LazyMultilinearInterpError <: Exception
    msg::String
end
function Base.showerror(io::IO, e::LazyMultilinearInterpError)
    print(io, "LazyMultilinearInterpError: ", e.msg)
end

"""
    lazy_multilinear_interpolation(params, nodes, grid; kwargs...)

This function is does multidimensional linear interpolation on a grid where the first two dimensions
represent the values being interpolated.  In other words, it interpolates matrices. It is written to
minimize the number of reads from the grid, and access things in an efficient order.
It is much faster than Interpolations.linear_interpolation when applied to memory-mapped grids.
At present it is used only for model atmosphere and departure coefficient interpolation.

See also: [`interpolate_marcs`](@ref), [`Korg.CubicSplines.CubicSpline`](@ref)
"""
function lazy_multilinear_interpolation(params, nodes, grid;
                                        param_names=["param $i" for i in 1:length(params)],
                                        perturb_at_grid_values=false)
    if perturb_at_grid_values
        # add small offset to each parameter which is exactly at grid value
        # this prevents the derivatives from being exactly zero
        on_grid_mask = in.(params, nodes)
        params[on_grid_mask] .= nextfloat.(params[on_grid_mask])

        # take care of the case where the parameter is at the last grid value
        too_high_mask = params .> last.(nodes)
        params[too_high_mask] .= prevfloat.(params[too_high_mask])
        params[too_high_mask] .= prevfloat.(params[too_high_mask])
    end

    upper_vertex = map(zip(params, param_names, nodes)) do (p, p_name, p_nodes)
        if !(p_nodes[1] <= p <= p_nodes[end])
            msg = "Can't interpolate grid.  $(p_name) is out of bounds. ($(p) ∉ [$(first(p_nodes)), $(last(p_nodes))])." *
                  " If you got the message calling `interpolate_marcs`, consider setting `clamp_abundances=true`."
            throw(LazyMultilinearInterpError(msg))
        end
        findfirst(p .<= p_nodes)
    end
    isexact = params .== getindex.(nodes, upper_vertex) #which params are on grid points?

    # allocate 2^n cube for each quantity
    dims = Tuple(2 for _ in upper_vertex) #dimensions of 2^n hypercube
    structure_type = typeof(promote(params...)[1])
    structure = Array{structure_type}(undef, (size(grid)[1:2]..., dims...))

    #put bounding atmospheres in 2^n cube
    for I in CartesianIndices(dims)
        local_inds = collect(Tuple(I))
        atm_inds = copy(local_inds)
        atm_inds[isexact] .= 2 #use the "upper bound" as "lower bound" if the param is on a grid point
        atm_inds .+= upper_vertex .- 2

        @views structure[:, :, local_inds...] .= grid[:, :, atm_inds...]
    end

    for i in eachindex(params) #loop over interpolation parameters
        isexact[i] && continue #no need to do anything for exact params

        # the bounding values of the parameter you are currently interpolating over 
        p1 = nodes[i][upper_vertex[i]-1]
        p2 = nodes[i][upper_vertex[i]]

        # inds1 and inds2 are the expressions for the slices through the as-of-yet 
        # uninterpolated quantities (temp, logPg, etc for atm itp) for each node value 
        # of the quantity being interpolated
        # inds1 = (1, 1, 1, ..., 1, 1, :, :, ...)
        # inds2 = (1, 1, 1, ..., 1, 2, :, :, ...)
        inds1 = vcat([1 for _ in 1:i-1], 1, [Colon() for _ in i+1:length(params)])
        inds2 = vcat([1 for _ in 1:i-1], 2, [Colon() for _ in i+1:length(params)])

        x = (params[i] - p1) / (p2 - p1) #maybe try using masseron alpha later
        @views structure[:, :, inds1...] = (1 - x) * structure[:, :, inds1...] +
                                           x * structure[:, :, inds2...]
    end
    structure[:, :, ones(Int, length(params))...]
end
