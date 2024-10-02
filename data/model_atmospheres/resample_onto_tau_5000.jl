# resample the built-it grid of marcs model atmospheres for cool dwarfs onto a grid of consistent 
# τ_5000 values. This was shown to improve accuracy tremendously in the Korg II paper.

# this script was used to generate the file, but it isn't part of the Korg module.  The file is 
# shipped via the artifact system alongside the big atmosphere archive.

using Korg, HDF5, ProgressMeter, Interpolations

function resample_grid(grid, new_log_taus; tau_index=4)
    dims = size(grid)[3:end] # number of m_h, alpha_m, c_m points
    resampled_grid = Array{Float32}(undef, (length(new_log_taus), 5, dims...))

    @showprogress for I in CartesianIndices(dims)
        for quant_ind in 1:5
            if quant_ind != tau_index
                itp = linear_interpolation(log10.(grid[:, tau_index, I]), grid[:, quant_ind, I];
                                           extrapolation_bc=Interpolations.Line())
                resampled_grid[:, quant_ind, I] = itp.(new_log_taus)
            end
        end
        resampled_grid[:, tau_index, I] .= 10 .^ new_log_taus

        # replace extrapolated values with nans.
        # layers filled with NaNs will be dropped when constructing the atmosphere object
        #taumin, taumax = grid[begin, tau_index, I], grid[end, tau_index, I]
        #resampled_grid[.! (log10(taumin) .< new_log_taus .< log10(taumax)), :, I] .= NaN
    end
    resampled_grid
end

nodes, grid = copy.(Korg._sdss_marcs_atmospheres);

# cool dwarfs only
Teff_mask = (nodes[1] .<= 4000) .& (nodes[1] .!= 3750)
logg_mask = nodes[2] .>= 3.5
nodes[1] = nodes[1][Teff_mask]
nodes[2] = nodes[2][logg_mask]
grid = grid[:, :, Teff_mask, logg_mask, :, :, :]

new_log_taus = -6:0.1:2

resampled_grid = resample_grid(grid, new_log_taus)

# save to disk
h5open("resampled_cool_dwarf_atmospheres.h5", "w") do f
    write_attribute(f, "version", "2024-01-04")
    f["grid"] = resampled_grid
    f["grid_parameter_names"] = ["Teff", "logg", "metallicity", "alpha", "carbon"]
    for i in 1:5
        f["grid_values/$i"] = nodes[i]
    end
end
