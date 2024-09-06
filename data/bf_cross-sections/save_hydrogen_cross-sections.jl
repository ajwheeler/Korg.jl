# run with
# julia --project save_hydrogen_cross-sections.jl path/to/Korg_data
using HDF5
include("parse_NORAD_H_cross_sections.jl")

DATA_DIR = ARGS[1]
OUT_FILE = "individual_H_cross-sections.h5"
cross_sections = parse_NORAD_H_cross_sections(joinpath(DATA_DIR,
                                                       "bf_cross_sections/NORAD/px.h1.1-40.txt"))

# total cross-sections for all orbitals for each value of n
ns = 1:40
Es = []
sigmas = []
for n in ns
    n_sigmas = filter(cross_sections) do xc
        first(xc) == n
    end
    _, _, _, E, _ = n_sigmas[1]
    total_n_sigma = zeros(length(E))
    for (_, L, Ebind, _, σs) in n_sigmas
        total_n_sigma .+= σs * (2 * (2L + 1))
    end
    push!(Es, E)
    push!(sigmas, total_n_sigma)
end

h5open(OUT_FILE, "w") do f
    f["n"] = collect(ns)
    f["E"] = hcat(Es...)
    f["sigma"] = hcat(sigmas...)
end
