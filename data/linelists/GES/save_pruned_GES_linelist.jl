using Korg

linelist = Korg.get_GES_linelist()

all_pruned = map([(6000, 4), (3000, 5), (4000, 3.5)]) do (Teff, logg)
    atm = interpolate_marcs(Teff, logg, A_X)
    Korg.prune_linelist(atm, linelist, A_X, wls; threshold=1e-4, sort_by_EW=false)
end

A_X = format_A_X()
wls = Korg.Wavelengths(linelist[1].wl, linelist[end].wl)

all_pruned = map([(6000, 4), (3000, 5), (4000, 3.5)]) do (Teff, logg)
    atm = interpolate_marcs(Teff, logg, A_X)
    Korg.prune_linelist(atm, linelist, A_X, wls; threshold=1e-4, sort_by_EW=false)
end

pruned = sort(unique([all_pruned...;]); by=l -> l.wl)
println(length(pruned) / length(linelist), " of the GES linelist was retained.")
