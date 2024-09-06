# Overview
#     Standalone script designed to produce plots to compare the calculation of the absorption
#     coefficient, using data from the Opacity Project, against absorption coefficients taken
#     from other sources. Currently, this just compares the absorption coefficients computed for
#     Hydrogenic atoms.

using Plots, Korg
include("../test/absorption_comparison_funcs.jl")

function plot_hydrogenic_comparison(spec, T=7800.0, ndens_species=3.0e16; p=nothing)
    λ_vals = 500:1.0:100_000
    hydrogenic_α_OP = OP_compare.calc_hydrogenic_bf_absorption_coef(λ_vals, T, ndens_species, spec;
                                                                    use_OP_data=true)
    hydrogenic_α_dflt = OP_compare.calc_hydrogenic_bf_absorption_coef(λ_vals, T, ndens_species,
                                                                      spec;
                                                                      use_OP_data=false)
    _pos_value_extrema(arr) = extrema(filter(x -> x >= 0, arr))
    ylims = extrema((_pos_value_extrema(hydrogenic_α_OP)...,
                     _pos_value_extrema(hydrogenic_α_dflt)...))
    display(ylims)
    if isnothing(p)
        p = plot()
    end
    plot!(p, λ_vals, [hydrogenic_α_OP, hydrogenic_α_dflt];
          labels=["Opacity Project" "Korg default (uses MHD)"],
          yscale=:log10, xscale=:log10, ylims=ylims,
          ylabel="absorption coef, α (cm⁻¹)", xlabel="λ (Ångstroms)",
          title="$spec Comparison")
    if spec == Korg.species"H I"
        vline!(p, [911.75, 3646, 8204, 14580, 22790, 32820]; label="known binding energies")
    end

    title!(p, "Opacity Project Comparison: $spec\nT = $T K, n($spec) = $ndens_species cm⁻³")
    return p
end

function main()
    T = 7800.0 # kind of arbitrary
    ndens_species = 3.0e16 # kind of arbitrary
    for spec in Korg.Species.(["H I"]) #He II no longer has a non-OP implementation.
        plot_fname = joinpath(@__DIR__, string(spec, " comp.png"))
        println("\nCreating Comparison plot at T = $T K, n($spec)= $ndens_species cm⁻³ for $spec")
        p = plot_hydrogenic_comparison(spec, T, ndens_species)
        println("Writing file to $plot_fname")
        savefig(p, plot_fname)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
