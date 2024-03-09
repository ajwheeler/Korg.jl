"""
TODO
"""
function weedout(atm, linelist, A_X, wls...; 
                 first_pass_threshold=1e-3, threshold=0.05, synthesis_kwargs...)

    # linelist will be sorted after call to synthesize
    sol = synthesize(atm, linelist, A_X, wls...; synthesis_kwargs...)
    wl_ranges = construct_wavelength_ranges(wls...)

    strong_lines = Line[]

    temps = get_temps(atm)
    β =  1 ./ (kboltz_eV * temps)
    n_div_Z = map(unique([l.species for l in linelist])) do spec
        spec => @. sol.number_densities[spec] / Korg.default_partition_funcs[spec](log(temps))
    end |> Dict

    λ_ind = 1
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

        E_upper = line.E_lower + c_cgs * hplanck_eV / line.wl 
        levels_factor = @. exp(-β*line.E_lower) - exp(-β*E_upper)
        # line center amplitude if it were a 1 Å tophat
        α_λ_line_center = @. 1e8 * 10.0^line.log_gf*sigma_line(line.wl)*levels_factor*n_div_Z[line.species]

        if any(α_λ_line_center .> first_pass_threshold .* sol.alpha[:, λ_ind])
            push!(strong_lines, line)
        end
    end

    really_strong_lines = Line[]
    @showprogress "checking $(length(strong_lines)) lines with synthesis" for line in strong_lines
        line_center = line.wl*1e8
        sol = synthesize(atm, [line], A_X, line_center, line_center; hydrogen_lines=false, synthesis_kwargs...)
        if 1 .- minimum(sol.flux ./ sol.cntm) > threshold
            push!(really_strong_lines, line)
        end
    end
    really_strong_lines
end