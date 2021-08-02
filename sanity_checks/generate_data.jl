using Korg, HDF5

#get output desination from command line
outfile = ARGS[1]

atm = read_model_atmosphere("../../Korg/test/data/sun.krz")
extralines = read_linelist("data/linelists/threshold01.vald")

#generate a rectified solar spectrum in a few key regions
regions = [
    #name, lower wl, upper wl, resolving power R, filenames
    ("Ca_H_K", 3930, 3970, 494000.0, "data/linelists/Ca-H-K" .* [".vald", "-2.vald"]),
    ("Mg_triplet", 5160, 5190, 432000.0, "data/linelists/Mg-triplet" .* [".vald", "-2.vald"]), 
    ("H_alpha", 6540, 6585, 698000.0, "data/linelists/Halpha" .* [".vald", "-2.vald", "-3.vald"])
   ]
for (name, lb, ub, R, fnames) in regions
    linelist = vcat(extralines, read_linelist.(fnames)...)
    sort!(linelist, by=line->line.wl)

    wls = lb-50:0.01:ub+50
    sol = synthesize(atm, linelist, wls; vmic=1.5)
    rF = Korg.rectify(Korg.constant_R_LSF(sol.flux, wls, R), wls)
    h5write(outfile, "spectra/sun/"*name, rF)
    h5writeattr(outfile, "spectra/sun/"*name, Dict("lb"=>lb, "ub"=>ub, "step"=>step(wls), "R"=>R))
end
