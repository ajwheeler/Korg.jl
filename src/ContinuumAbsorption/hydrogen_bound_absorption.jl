using Interpolations: LinearInterpolation, Line
using HDF5

# load hydrogen bf cross sections
const _H_cross_sections = let
    h5open(joinpath(_data_dir, "bf_cross-sections", 
                                         "individual_H_cross-sections.h5")) do f
        sigmas = map(eachcol(read(f["E"])), eachcol(read(f["sigma"]))) do Es, σs
            LinearInterpolation(Es, σs, extrapolation_bc=Line())
        end
        # use the cross sections for the first 6 energy levels only.
        # the binding energy for n=7 corresponds to ~45,000 Å 
        collect(zip(read(f["n"]), sigmas))
    end
end
