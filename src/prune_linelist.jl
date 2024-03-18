using Statistics: mean

"""
    prune_linelist(atm, linelist, A_X, wls...; threshold=1.0, sort=true, synthesis_kwargs...)

Return the vector containing the strongest lines in  `linelist`, (optionally) sorted by approximate 
equivalent width.

# Arguments
- `atm`: the atmosphere model
- `linelist`: the linelist (a vector of `Line` objects)
- `A_X`: the abundance of each element (see [`format_A_X`](@ref))
- `wls...`: the wavelength ranges to synthesize over.  These are specified the same way as the 
  `wls` for [`synthesize`](@ref).

# Keyword Arguments
- `threshold=0.1`: The threshold ratio in the line center absorption to the continuum absorption 
   computed at the photosphere for a line to be included in the returned list.
   `0.1` is a reasonable default for getting a sense of what might be measurable in a high-res, 
   high-quality spectrum, but it should not be used to create a linelist for synthesis.  
- `sort_by_EW=true`: If `true`, the returned linelist will be sorted by approximate equivalent 
   width.  If `false`, the linelist will be in wavelength order. Leaving the list in wavelength 
   order is much faster, but sorting by strength is useful for visualizing the strongest lines.
All other kwargs are passed to internal calls to [`synthesize`](@ref).

!!! caution
    While this function can be used to prune a linelist for synthesis, the default behavior too 
    aggressive for this purpose.  Set a much lower threshold (e.g. `threshold=1e-4`) and use 
    `sort=false` if you are pruning the linelist to speedup synthesis.  Note that Korg will 
    dynamically choose which lines to include even if you use a large linelist (see 
    the `line_cutoff_threshold` keyword argument to [`synthesize`](@ref)).

See also [`merge_close_lines`](@ref) if you are using this for plotting.
"""
function prune_linelist(atm, linelist, A_X, wls...; 
                        threshold=0.1, sort_by_EW=true, synthesis_kwargs...)
    # linelist will be sorted after call to synthesize
    sol = synthesize(atm, linelist, A_X, wls...; synthesis_kwargs...)
    cntm_sol = synthesize(atm, [], A_X, wls...; synthesis_kwargs...) 
    wl_ranges = construct_wavelength_ranges(wls...)

    # get the atmosphere layer where τ ≈ 1 for each wavelength
    approximate_τ = cumsum(sol.alpha[1:end-1, :] .* -diff(get_zs(atm)), dims=1)
    photosphere_indices = map(eachcol(approximate_τ)) do τ
        findfirst(τ .> 1.0)
    end

    # precompute T, thermo β, n/Z, and Doppler widths at each level
    temps = get_temps(atm)
    β =  1 ./ (kboltz_eV * temps)
    linelist_species = unique([l.species for l in linelist])
    n_div_Z = map(linelist_species) do spec
        spec => @. sol.number_densities[spec] / Korg.default_partition_funcs[spec](log(temps))
    end |> Dict
    doppler_widths = map(linelist_species) do spec
        spec => doppler_width.(sol.wavelengths[1]*1e-8, get_temps(atm), get_mass(spec), 0.0) * 1e8
    end |> Dict

    λ_ind = 1
    strong_lines = Line[]
    for line in linelist
        line_center = line.wl * 1e8
        if !any(λs[begin] < line_center < λs[end] for λs in wl_ranges)
            continue
        end

        # move λ_ind to the wavelength in the synthesis grid closest to the line center
        while (λ_ind < length(sol.wavelengths)) && 
              abs(sol.wavelengths[λ_ind] - line_center) > abs(sol.wavelengths[λ_ind + 1] - line_center)
            λ_ind += 1
        end
        phot_ind = photosphere_indices[λ_ind]

        E_upper = line.E_lower + c_cgs * hplanck_eV / line.wl 
        levels_factor = exp(-β[phot_ind]*line.E_lower) - exp(-β[phot_ind]*E_upper)
        # line center amplitude if it were a 1 Å tophat

        σ = doppler_widths[line.species][phot_ind]
        n_Z = n_div_Z[line.species][phot_ind]
        α_λ_line_center = 1e8 * 10.0^line.log_gf*sigma_line(line.wl)*levels_factor*n_Z / σ

        if α_λ_line_center > threshold * cntm_sol.alpha[phot_ind, λ_ind]
            push!(strong_lines, line)
        end
    end

    # sort lines by approximate EW or leave them in wavelength order
    if sort_by_EW
        approx_EWs = @showprogress "measuring $(length(strong_lines)) lines" map(strong_lines) do line
            line_center = line.wl*1e8
            sol = synthesize(atm, [line], A_X, line_center - 2.0, line_center + 2.0; 
                             hydrogen_lines=false, use_chemical_equilibrium_from=sol, synthesis_kwargs...)
            sum(1 .- sol.flux ./ sol.cntm) # units don't matter
        end
        strong_lines[sortperm(approx_EWs, rev=true)]
    else
        strong_lines
    end
end

"""
    merge_close_lines(linelist; merge_distance=0.2)

Produce a list of the species and wavelengths of the lines in `linelist`, merging lines of the same 
species that are within `merge_distance` (default: 0.2 Å).  This is useful for labeling lines in a 
plot after running [`prune_linelist`](@ref).

# Arguments
- `linelist`: the linelist (a vector of `Line` objects)

# Keyword Arguments
- `merge_distance=0.2`: The maximum distance in Å between lines of the same species to be merged 
   into a single entry

# Returns
A vector of tuples `(wl, species)` where `wl` is the wavelength of the line in Å and `species` is
a string identifying the species of the line.
"""
function merge_close_lines(lines; merge_distance=0.2)
    lines = sort(lines, by=l->l.wl)
    
    all_species = unique(l.species for l in lines)
    species_line_lists = map(all_species) do species
        species_lines = filter(l -> l.species == species, lines)

        merged_lines = []
        to_merge = [species_lines[1]]

        for line in species_lines[2:end]
            if (line.wl - to_merge[end].wl) * 1e8 < merge_distance
                push!(to_merge, line)
            else
                push!(merged_lines, (mean(l.wl*1e8 for l in to_merge), string(species)))
                to_merge = [line]
            end
        end
        push!(merged_lines, (mean(l.wl*1e8 for l in to_merge), string(species)))
        merged_lines
    end

    sort(vcat(species_line_lists...))
end