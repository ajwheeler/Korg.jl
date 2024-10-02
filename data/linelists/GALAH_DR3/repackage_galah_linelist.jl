# this script repackages the GALAH DR3 linelist into an HDF5 file containing only the columns used 
# by Korg

using FITSIO, HDF5, DataFrames, Korg

# this should point to the base dir of the Korg_data repo
# https://github.com/ajwheeler/Korg_data
Korg_data_dir = "../../../../Korg_data"

f = FITS(joinpath(Korg_data_dir, "linelists", "GALAH_dr3", "galah_master_v5.2.fits"))
galah_lines = DataFrame(f[2]);
galah_lines.wl = Korg.air_to_vacuum.(galah_lines.LAMBDA)
galah_lines.species = map(galah_lines.NAME, galah_lines.ION) do name, ion
    join([name; " "; string(ion)])
end

galah_lines.formula = map(galah_lines.NAME) do name
    map(name) do el
        if el == " "
            0x0
        else
            Korg.atomic_numbers[el]
        end
    end
end

# calculate isotopic corrections to log gfs
Δlog_gf = map(eachrow(galah_lines)) do line
    f = 1.0
    for (el, iso) in zip(line.NAME, line.ISOTOPE)
        if el in Korg.atomic_symbols && iso != 0
            Z = Korg.atomic_numbers[el]
            f *= Korg.isotopic_abundances[Z][iso]
        end
    end
    log10(f)
end
galah_lines.LOG_GF .+= Δlog_gf

h5open("galah_dr3_linelist.h5", "w") do f
    f["formula", compress=9] = hcat(galah_lines.formula...)
    f["ionization", compress=9] = galah_lines.ION
    f["wl", compress=9] = Float32.(galah_lines.wl)
    f["log_gf", compress=9] = Float32.(galah_lines.LOG_GF)
    f["E_lo", compress=9] = Float32.(galah_lines.E_LOW)
    f["gamma_rad", compress=9] = Float32.(galah_lines.RAD_DAMP)
    f["gamma_stark", compress=9] = Float32.(galah_lines.STARK_DAMP)
    f["vdW", compress=9] = Float32.(galah_lines.VDW_DAMP)
end
