using Korg, HDF5, FITSIO

# this should point to the base dir of the Korg_data repo
# https://github.com/ajwheeler/Korg_data
Korg_data_dir = joinpath(@__DIR__, "../../../../Korg_data")
linelist_dir = joinpath(Korg_data_dir, "linelists/GES")

atomic_lines = FITS(joinpath(linelist_dir, "J_A+A_645_A106_geslines.dat.gz.fits")) do f
    map(read(f[2], "Element"),
        read(f[2], "Ion"),
        read(f[2], "Isotope"),
        read(f[2], "lambda"),
        read(f[2], "loggf"),
        read(f[2], "Elow"),
        read(f[2], "Rad-damp"),
        read(f[2], "Sta-damp"),
        read(f[2], "Vdw-damp")) do el, ion, isotope, wl, loggf, Elow, gamma_rad, gamma_stark,
                                   gamma_vdw
        spec = Korg.Species(Korg.Formula(el), ion - 1)
        Z = Korg.atomic_numbers[el]

        Δloggf = if Z in keys(Korg.isotopic_abundances) &&
                    isotope in keys(Korg.isotopic_abundances[Z])
            log10(Korg.isotopic_abundances[Z][isotope])
        else
            0
        end

        (spec, wl, Float32(loggf + Δloggf), Elow, gamma_rad, gamma_stark, gamma_vdw)
    end
end
# Remove hydrogren lines.  There are no triply-ionized species in the GES linelist to worry about.
filter!(atomic_lines) do line
    Korg.get_atoms(line[1]) != [1]
end

mol_lines = FITS(joinpath(linelist_dir, "J_A+A_645_A106_gesmol.dat.gz.fits")) do f
    map(read(f[2], "Element1"),
        read(f[2], "Element2"),
        read(f[2], "Isotope1"),
        read(f[2], "Isotope2"),
        read(f[2], "lambda"),
        read(f[2], "loggf"),
        read(f[2], "Elow"),
        read(f[2], "Rad-damp")) do el1, el2, isotope1, isotope2, wl, loggf, Elow, gamma_rad
        spec = Korg.Species(el1 * el2)
        Z1 = Korg.atomic_numbers[el1]
        Z2 = Korg.atomic_numbers[el2]

        Δloggf = map([(Z1, isotope1), (Z2, isotope2)]) do (Z, isotope)
            if Z in keys(Korg.isotopic_abundances) &&
               isotope in keys(Korg.isotopic_abundances[Z])
                log10(Korg.isotopic_abundances[Z][isotope])
            else
                0.0f0
            end
        end |> sum

        (spec, wl, Float32(loggf + Δloggf), Elow, gamma_rad, 0.0f0, 0e0)
    end
end

linelist = sort([atomic_lines; mol_lines]; by=l -> l[2])

h5open("Heiter_et_al_2021.h5", "w") do f
    f["species"] = string.(first.(linelist))
    f["wl", compress=3] = [l[2] for l in linelist]
    f["log_gf", compress=3] = [l[3] for l in linelist]
    f["E_lower", compress=3] = [l[4] for l in linelist]
    f["gamma_rad", compress=3] = [l[5] for l in linelist]
    f["gamma_stark", compress=3] = [l[6] for l in linelist]
    f["vdW"] = [l[7] for l in linelist]
end
