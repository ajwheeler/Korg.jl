using HDF5, Korg

# this should point to the base dir of the Korg_data repo
# https://github.com/ajwheeler/Korg_data
Korg_data_dir = "../../../../Korg_data"

filepath = joinpath(Korg_data_dir, "linelists", "APOGEE_DR17", "turbospec.h2o_POKAZATEL-BC9.5.dat")
water = read_linelist(filepath, format="turbospectrum")

h5open("pokazatel_water_lines.h5", "w") do f
    f["wl", compress=3] = [Float32(l.wl) for l in water]
    f["log_gf", compress=3] = [Float32(l.log_gf) for l in water]
    f["E_lower", compress=3] = [Float32(l.E_lower) for l in water]
    f["gamma_rad", compress=3] = [Float32(l.gamma_rad) for l in water]
end