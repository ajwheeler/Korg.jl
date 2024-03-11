using Statistics: mean

"""
TODO
"""
function weedout(atm, linelist, A_X, wls...; 
                 first_pass_threshold=1e-3, threshold=0.1, synthesis_kwargs...)

    # linelist will be sorted after call to synthesize
    sol = synthesize(atm, linelist, A_X, wls...; synthesis_kwargs...)
    cntm_sol = synthesize(atm, [], A_X, wls...; synthesis_kwargs...) 
    wl_ranges = construct_wavelength_ranges(wls...)

    approximate_τ = cumsum(sol.alpha[1:end-1, :] .* -diff(get_zs(atm)), dims=1)
    photosphere_indices = map(eachcol(approximate_τ)) do τ
        findfirst(τ .> 1.0)
    end

    σs = doppler_width.(sol.wavelengths[1]*1e-8, get_temps(atm), atomic_masses[1], 0.0) * 1e8

    strong_lines = Line[]

    temps = get_temps(atm)
    β =  1 ./ (kboltz_eV * temps)
    n_div_Z = map(unique([l.species for l in linelist])) do spec
        spec => @. sol.number_densities[spec] / Korg.default_partition_funcs[spec](log(temps))
    end |> Dict

    λ_ind = 1
    for line in linelist
        isna = (line.species == species"Na") && (5891 < line.wl*1e8 < 5892) && (line.log_gf > -1)
        line_center = line.wl * 1e8
        if !any(λs[begin] < line_center < λs[end] for λs in wl_ranges)
            continue
        end
        isna && println(line_center)

        # move λ_ind to the wavelength in the synthesis grid closest to the line center
        while (λ_ind < length(sol.wavelengths)) && 
              abs(sol.wavelengths[λ_ind] - line_center) > abs(sol.wavelengths[λ_ind + 1] - line_center)
            λ_ind += 1
        end
        phot_ind = photosphere_indices[λ_ind]


        isna && println(phot_ind)
        isna && println(approximate_τ[:, λ_ind])
        isna && println(findfirst(approximate_τ[:, λ_ind] .> 1))
        isna && println()

        E_upper = line.E_lower + c_cgs * hplanck_eV / line.wl 
        levels_factor = exp(-β[phot_ind]*line.E_lower) - exp(-β[phot_ind]*E_upper)
        # line center amplitude if it were a 1 Å tophat
        α_λ_line_center = 1e8 * 10.0^line.log_gf*sigma_line(line.wl)*levels_factor*n_div_Z[line.species][phot_ind] / σs[phot_ind]

        isna && println(α_λ_line_center, " ", cntm_sol.alpha[phot_ind, λ_ind])
        isna && println()

        if α_λ_line_center > first_pass_threshold * cntm_sol.alpha[phot_ind, λ_ind]
            push!(strong_lines, line)
        end
    end

    really_strong_lines = Line[]
    @showprogress "checking $(length(strong_lines)) lines with synthesis" for line in strong_lines
        #line_center = line.wl*1e8
        #sol = synthesize(atm, [line], A_X, line_center, line_center; hydrogen_lines=false, synthesis_kwargs...)
        #if 1 .- minimum(sol.flux ./ sol.cntm) > threshold
            push!(really_strong_lines, line)
        #end
    end
    really_strong_lines
end

"""
TODO
"""
function merge_close_lines(lines; merge_distance=1.0)
    @assert issorted(lines, by=l->l.wl)
    
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