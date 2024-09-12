import Pkg;
using Dates;
t0 = now();
t_then = t0;
using InteractiveUtils;
versioninfo();
Pkg.activate("/n/home12/saydjari/finksage/julia_env/korg_1p10/");
Pkg.instantiate();
Pkg.precompile(); ## FOR ZACK: You will need to change this path to your own envionrment
t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Package activation took $dt");
t_then = t_now;
flush(stdout);

using Distributed, SlurmClusterManager, Suppressor, DataFrames
addprocs(SlurmManager())

t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Worker allocation took $dt");
t_then = t_now;
flush(stdout);

@everywhere begin
    using ProgressMeter, Korg, Random, DelimitedFiles, Distributions, HDF5, ParallelDataTransfer,
          ThreadPinning
    using LinearAlgebra, BLISBLAS, FITSIO

    runName = "KommenceKorg"
    # B 3600 to 5800
    # R 5760 to 7620
    # Z 7520 to 9824
    wl_lo, wl_hi = 3000, 9000 # the wavelength grid is 3000 : 0.01 : 9000
    # Do we we want to pad this? Or force the resolution to be lower? Convolution either way.
    desilines = Korg.get_VALD_solar_linelist() ## FOR ZACK: We will want to replace this at some point with a better linelist
    # LL is 3000 to 9000

    # get the parameter ranges within with the SDSS atmosphere grid exists
    grid_vals = Korg.get_atmosphere_archive()[1]
    atm_lbs = first.(grid_vals)
    atm_ubs = last.(grid_vals)

    function query_spectra(intup)
        indx, Teff, logg, A_Xi, vmic = intup
        A_X = collect(A_Xi)

        try
            atm = Korg.interpolate_marcs(Teff, logg, A_X;)
            tout = synthesize(atm, desilines, A_X, wl_lo, wl_hi; vmic=vmic, hydrogen_lines=true,
                              hydrogen_line_window_size=300,
                              electron_number_density_warn_threshold=1e10)

            dirindx = lpad(indx ÷ 1000, 3, "0")
            specindx = lpad(indx, 6, "0")
            dir = "$(runName)/$(dirindx)"
            mkpath(dir)
            writedlm(joinpath(dir, "spectrum_$specindx.dat"), [tout.flux (tout.flux) ./ tout.cntm])
            return true
        catch
            return false
        end
        GC.gc()
    end
end
Pkg.status("Korg")
t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Korg atm loading took $dt");
t_then = t_now;
flush(stdout);
println(BLAS.get_config());
flush(stdout);

rng = MersenneTwister(368)

## Uniform Prior Section
Nspec = 5_000
# A_X[1] is definitionally fixed to 12 and is for H
A_X = repeat(Korg.grevesse_2007_solar_abundances'; inner=(Nspec, 1))
A_X[:, 2:end] .+= 0.01 * randn(rng, Nspec, Korg.MAX_ATOMIC_NUMBER - 1); #0.01 dex sigma

m_h = rand(rng, Uniform(atm_lbs[3], atm_ubs[3]), Nspec)
A_X[:, 3:Korg.MAX_ATOMIC_NUMBER] .= A_X[:, 3:Korg.MAX_ATOMIC_NUMBER] .+ m_h

a_m = rand(rng, Uniform(atm_lbs[4], atm_ubs[4]), Nspec)
for alpha_el in [8, 12, 14, 20, 22] # O Mg Si Ca Ti
    A_X[:, alpha_el] = A_X[:, alpha_el] .+ a_m
end
c_m = zeros(Nspec)

Teff = rand(rng, Uniform(atm_lbs[1], atm_ubs[1]), Nspec)
logg = rand(rng, Uniform(atm_lbs[2], atm_ubs[2]), Nspec)
vmic_min = 0
vmic_max = 5
vmic = rand(rng, Uniform(vmic_min, vmic_max), Nspec);

println("Teff Uniform Prior: $(atm_lbs[1]), $(atm_ubs[1])")
println("logg Uniform Prior: $(atm_lbs[2]), $(atm_ubs[2])")
println("m_h Uniform Prior: $(atm_lbs[3]), $(atm_ubs[3])")
println("a_m Uniform Prior: $(atm_lbs[4]), $(atm_ubs[4])")
println("vmic Uniform Prior: $vmic_min, $vmic_max")

mkpath(runName)
h5write("$(runName)/korg_grid_params.h5", "Teff", Teff)
h5write("$(runName)/korg_grid_params.h5", "logg", logg)
h5write("$(runName)/korg_grid_params.h5", "m_h", m_h)
h5write("$(runName)/korg_grid_params.h5", "a_m", a_m)
h5write("$(runName)/korg_grid_params.h5", "c_m", c_m)
h5write("$(runName)/korg_grid_params.h5", "vmic", vmic)
h5write("$(runName)/korg_grid_params.h5", "abundances", A_X)

itarg = Iterators.zip(1:Nspec, Teff, logg, eachrow(A_X), vmic);

t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Preparing for and launching pmap $dt");
t_then = t_now;
flush(stdout);
pout = @showprogress pmap(query_spectra, itarg);
h5write("$(runName)/korg_grid_params.h5", "converged_flag", convert.(Int, pout))
rmprocs(workers())