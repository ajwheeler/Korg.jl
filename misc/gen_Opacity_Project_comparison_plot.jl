# Overview
#     Standalone script designed to produce plots to compare the calculation of the absorption
#     coefficient, using data from the Opacity Project, against absorption coefficients taken
#     from other sources. Currently, this just compares the absorption coefficients computed for
#     Hydrogenic atoms.

using Plots, Korg
include("../test/opacity_comparison_funcs.jl")

function _pos_value_extrema(arr)
    min_val = floatmax(eltype(arr))
    max_val = 0.0

    for elem in arr
        if !(elem > 0)
            continue
        end
        min_val = min(min_val, elem)
        max_val = max(max_val, elem)
    end
    (min_val,max_val)
end

function _plot_x_intervals!(p, x_intervals)
    y_lim = ylims(p)
    for i in 1:size(x_intervals)[1]
        cur_x = [x_intervals[i,1], x_intervals[i,2]]
        cur_y1 = [y_lim[1], y_lim[1]]
        cur_y2 = [y_lim[2], y_lim[2]]
        plot!(p, cur_x, cur_y1, fillrange = cur_y2, fillalpha = 0.10,
              ylim = y_lim, fillcolor = :grey,
              label = "")
    end
end

function plot_hydrogenic_comparison(species_name, T = 7800.0, ndens_species = 3.0e16; p = nothing)
    λ_vals = OP_compare._dflt_λ_vals(species_name)
    hydrogenic_α_OP = OP_compare.calc_hydrogenic_bf_absorption_coef(λ_vals,  T, ndens_species,
                                                                    species_name;
                                                                    use_OP_data = true)
    hydrogenic_α_dflt = OP_compare.calc_hydrogenic_bf_absorption_coef(λ_vals,  T, ndens_species,
                                                                      species_name;
                                                                      use_OP_data = false)
    ylims = extrema((_pos_value_extrema(hydrogenic_α_OP)...,
                     _pos_value_extrema(hydrogenic_α_dflt)...))
    if isnothing(p)
        p = plot()
    end
    p = plot!(p, λ_vals, [hydrogenic_α_OP, hydrogenic_α_dflt],
              labels = ["Opacity Project" "Reference (Hydrogenic)"],
              yscale = :log10, xscale = :log10, ylims = ylims,
              ylabel = "absorption coef, α (cm⁻¹)", xlabel = "λ (Ångstroms)",
              title = string(join(split(species_name,'_'), " "), " Comparison"))
    if species_name == "H_I"
        vline!(p, [911.75, 3646, 8204, 14580, 22790, 32820], label = "known binding energies")
    end
    λ_comp_intervals = OP_compare._λ_comp_intervals(species_name)
    _plot_x_intervals!(p, λ_comp_intervals)

    cleaned_name = join(split(species_name,'_'), " ")
    title_str = string("Opacity Project Comparison: $cleaned_name\n",
                       "T = $T K, n($cleaned_name) = $ndens_species cm⁻³")
    title!(p,title_str)

    return p
end


function main()
    T = 7800.0 # kind of arbitrary
    ndens_species = 3.0e16 # kind of arbitrary
    for species_name in ["H_I", "He_II"]
        plot_fname = joinpath(@__DIR__, string(species_name, "_comp.png"))
        println("\nCreating Comparison plot at T = $T K, ndens_species = $ndens_species cm⁻³ for ",
                species_name)
        p = plot_hydrogenic_comparison(species_name, T, ndens_species);
        println("Writing file to $plot_fname")
        savefig(p, plot_fname)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
