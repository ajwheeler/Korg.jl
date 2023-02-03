using HDF5
include("parse_NORAD_H_cross_sections.jl")

DATA_DIR = ARGS[1]
OUT_FILE = "individual_H_cross-sections.h5"
cross_sections = parse_NORAD_H_cross_sections(joinpath(DATA_DIR, "bf_cross_sections/NORAD/px.h1.1-40.txt"))

h5open(OUT_FILE, "w") do f
    f["n"] = first.(cross_sections)
    f["L"] = (t->t[2]).(cross_sections)
    f["E_bind"] = (t->t[3]).(cross_sections)
    f["E"] = hcat((t->t[4]).(cross_sections)...)
    f["sigma"] = hcat((t->t[4]).(cross_sections)...)
end