# Overview
#     Standalone script designed to produce a plot comparing the tabulated data from the
#     panels of Figure 8.5 of Gray's 2005 edition of "The Observation and Analysis of Stellar
#     Photospheres".
#
#     Data was extracted using the WebPlotDigitizer tool. We've recorded values with greater
#     precision than the values should actually have. The values are primarily meant to be used to
#     make sure that our implementation give results to the correct order of magniutde and
#     qualitatively has the expected wavelength dependences.
#
#     This isn't a particularly idiomatic approach, but it should get the job done
#
# Synopsis
#     `julia gen_opacity_plot.jl PANEL OUTPUT_NAME`
#
# Description
#     PANEL should be a, b, c, or d
#     OUTPUT_NAME the location where the opacity will be saved

using Plots, Korg
include("../test/absorption_comparison_funcs.jl")

function gen_plot(plot_fname, temperature, Pₑ, panel_dict)
    l = @layout [a; c{0.3h}]
    p = plot(; layout=l)
    for (i, (key, (orig_λ_vals, orig_data))) in enumerate(panel_dict)
        plot!(p[1], orig_λ_vals, orig_data; label=key,
              color=i)
        if haskey(Gray_opac_compare.Gray05_opacity_form_funcs, key)
            func, bounds, label_suffix = Gray_opac_compare.Gray05_opacity_form_funcs[key]
            w = Korg.ContinuumAbsorption.contained.(orig_λ_vals, Ref(bounds))
            λ_vals = orig_λ_vals[w]
            calculated_vals = func.(λ_vals, temperature, Pₑ)
            label = nothing #string("calculated ", label_suffix)
            plot!(p[1], λ_vals, calculated_vals * 1e26;
                  seriestype=:scatter, label=label, color=i)

            plot!(p[2], λ_vals, calculated_vals * 1e26 - orig_data[w];
                  color=i, label=nothing)
        end
    end
    ylabel!(p[1], "κ_ν ρ / (n(H I) Pₑ)\n[cm⁴ dyne⁻¹ per H atom]")
    ylabel!(p[2], "residual")
    xlabel!(p[2], "λ (Ångstroms)")

    savefig(p, plot_fname)
end

function main(args)
    # first, parse arguments
    if length(args) != 2
        msg = """
This script expects 2 arguments. The first argument should indicate the panel to compare data from.
The second argument should specify the file where the resulting figure should be saved.  The panels 
contain:
    a: H⁻ bf
    b: H⁻ ff, He⁻ ff, H⁻ bf, H₂⁺, H
    c: H⁻ ff, He⁻ ff, H⁻ bf, H I
    d: H I bf and ff
"""
        throw(ArgumentError(msg))
    end
    panel_name = args[1]
    plot_fname = args[2]
    if !(panel_name in ["a", "b", "c", "d"])
        throw(ArgumentError("Invalid panel name. The only valid names are a, b, c, or d"))
    end
    println("Plotting data for panel \"", panel_name, "\" and writing the result to ", plot_fname)

    # load the data extracted from the figure
    panel_prop, panel_dict = Gray_opac_compare.load_panel_data(panel_name)

    temperature = panel_prop["temperature"]
    Pₑ = 10^panel_prop["log10(electron_pressure)"]

    gen_plot(plot_fname, temperature, Pₑ, panel_dict)
end

main(ARGS)
