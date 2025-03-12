using Korg

linelist = Korg.get_GES_linelist()
A_X = format_A_X()
teff_loggs = [(6000, 4), (3000, 5), (4000, 2)]
atms = map(teff_loggs) do (Teff, logg)
    interpolate_marcs(Teff, logg, A_X)
end
wls = Korg.Wavelengths(linelist[1].wl, linelist[end].wl)

all_pruned = map(atms) do atm
    @time Korg.prune_linelist(atm, linelist, A_X, wls;
                              threshold=1e-6, sort_by_EW=false, max_distance=10)
end

pruned = sort(unique([all_pruned...;]); by=l -> l.wl)
println(length(pruned) / length(linelist), " of the GES linelist was retained.")

using PyPlot

map(atms, teff_loggs) do atm, (Teff, logg)
    figure(; figsize=(12, 5))
    @time spec1 = synthesize(atm, linelist, A_X, wls)
    @time spec = synthesize(atm, pruned, A_X, wls)
    title("Teff = $Teff, logg = $logg")
    plot(spec.wavelengths, (spec.flux - spec1.flux) ./ spec.cntm)
    ylabel("rect flux residuals")
    xlabel("Wavelength (Å)")
    savefig("GES_pruned_T_$(Teff)_logg_$(logg).png")
end