# This script is an example of how to use Korg to generate a grid of synthetic spectra using Slurm.
using Distributed, SlurmClusterManager, DataFrames, Distributions, HDF5, ProgressMeter

# This will automatically add a worker for each task.
addprocs(SlurmManager())
# if you you are using an project/environment (recommended), you can make sure it's activated on the
# workers like so:
# addprocs(SlurmManager(), exeflags = ["--project=."])) # or wherever your Project.toml is

@everywhere begin
    using Korg, Random, DelimitedFiles
    runName = "Korg-Slurm-Example"
end
println("Run name: $runName")
mkdir(runName) # do this now so that it errors fast if the directory already exists
@everywhere begin
    wl_lo, wl_hi = 5600, 5900 # wavelength range to synthesize
    linelist = Korg.get_GALAH_DR3_linelist()

    function generate_spectrum(spec_params)
        indx, Teff, logg, A_Xi = spec_params
        A_X = collect(A_Xi)
        atm = Korg.interpolate_marcs(Teff, logg, A_X)

        try
            # turn off the warning about electron number density because it will be thrown for
            # some unphysical Teff, logg combinations.
            sol = synthesize(atm, linelist, A_X, (wl_lo, wl_hi);
                             electron_number_density_warn_threshold=Inf)

            # write flux and rectified flux to file
            # place the files in directories of 1000 spectra each
            dirindx = lpad(indx ÷ 1000, 3, "0")
            specindx = lpad(indx, 6, "0")
            dir = "$(runName)/$(dirindx)"
            mkpath(dir)
            writedlm(joinpath(dir, "spectrum_$specindx.dat"), [sol.flux (sol.flux) ./ sol.cntm])

            return true
        catch e
            println("Failed to generate spectrum $indx")
            println(e)
            return false
        end
    end
end

# set of the the stellar parameters and abundances for each spectrum.
# Here, we will generate 100 spectra with random stellar parameters and abundances, but
# you may want a uniform grid, or something else.

nspec = 5000 # Number of spectra to generate
rng = MersenneTwister(368) # n.b. this is not guaranteed to be stable across Julia versions

# Here, we will choose a random Teff, logg, and [m/H] foreach spectrum, leaving everything else fixed.
# For real science, you will probably want to vary additional parameters (e.g. vmic, vsini,
# additional abundances). See the Korg documentation for how to do this.
# Note that many of these Teff/logg combinations will be unphysical.
Teff = rand(rng, Uniform(2800, 8000), nspec)
logg = rand(rng, Uniform(-0.5, 5.5), nspec)
M_H = rand(rng, Uniform(-2.5, 0.5), nspec)

# turn those metallicities into full abundance vectors
A_X = Matrix{Float64}(undef, Korg.MAX_ATOMIC_NUMBER, nspec)
for i in 1:nspec
    A_X[:, i] .= format_A_X(M_H[i])
end

h5write("$(runName)/korg_grid_params.h5", "Teff", Teff)
h5write("$(runName)/korg_grid_params.h5", "logg", logg)
h5write("$(runName)/korg_grid_params.h5", "M_H", M_H)
h5write("$(runName)/korg_grid_params.h5", "abundances", A_X)

# iterate over each combination of parameters, producing a spectrum
all_params = Iterators.zip(1:nspec, Teff, logg, eachcol(A_X));
success = @showprogress pmap(generate_spectrum, all_params);

h5write("$(runName)/korg_grid_params.h5", "converged_flag", convert.(Int, success))

rmprocs(workers()) # remove workers so slurm doesn't have to
