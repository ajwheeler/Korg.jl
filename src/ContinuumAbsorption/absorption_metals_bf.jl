metal_bf_cross_sections = let 
    cross_sections = Dict{Species, Any}()
    h5open(joinpath(_data_dir, "metal_bf_cross-sections", "bf_cross-sections.h5")) do f
        T_grid = read(f["logT_min"]) : read(f["logT_step"]) : read(f["logT_max"])
        ν_grid = read(f["nu_min"]) : read(f["nu_step"]) : read(f["nu_max"])
        for dataset in f["cross-sections"]
            #parse the dataset name to a species
            spec = Species(split(HDF5.name(dataset), "/")[end])
            σs = read(dataset)
            cross_sections[spec] = LinearInterpolation((ν_grid, T_grid), σs, extrapolation_bc=Flat())
        end
    end
    cross_sections
end

"""
    metal_bf_absorption!(α, νs, T, number_densities)

Adds to α the contributions of bf metal opacities.  Uses precomputed tables from 
[TOPBase](http://cdsweb.u-strasbg.fr/topbase/topbase.html) for Li, Be, B, C, N, O, F, Ne, Na, Mg, 
Al, Si, S, Ar, and Ca. Uses tables from 
[NORAD](https://www.astronomy.ohio-state.edu/nahar.1/nahar_radiativeatomicdata/index.html) for 
Fe.  For these elements, tables have been precomputed for the neutral and singly ionized species 
assuming and LTE distribution of energy levels.  See `Korg/data/metal_bf_cross-sections/` for the
scripts which generate the tables.

Cross sections were computed for 100 K < T < 100,000 K and frequencies corresponding to 
500 Å < λ < 30,000 Å.  Outside of either of those ranges, flat extrapolation is used (i.e. the 
extreme value is returned).
"""
function metal_bf_absorption!(α, νs, T, number_densities)
    for spec in keys(metal_bf_cross_sections)
        if spec in keys(number_densities)
            if spec in [species"H I", species"He I", species"H II"]
                continue #these are handled with other approximations
            end
            # cross-section interpolator
            σ_itp = metal_bf_cross_sections[spec]
            α .+= exp.(log(number_densities[spec]) .+ σ_itp.(νs, log10(T))) * 1e-18 #convert to cm^2 
        end
    end

end