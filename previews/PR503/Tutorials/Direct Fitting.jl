# [Download this page as a Jupyter notebook](./Direct Fitting.ipynb) #md
# You may have to execute these cells to see the plot outputs #nb
#
# In this example, we'll show how to use [`Korg.Fit.fit_spectrum`](@ref) to fit one of the spectra from
# [Griffith et al. 2022](https://ui.adsabs.harvard.edu/abs/2022arXiv221001821G). You will probably
# find it helpful to look at the [documentation for this
# function](https://ajwheeler.github.io/Korg.jl/stable/API/#Korg.Fit.fit_spectrum) as well. For
# fitting equivalent widths (rather than spectra directly) see [the documentation for
# `Korg.Fit.ews_to_abundances`](https://ajwheeler.github.io/Korg.jl/stable/API/#Korg.Fit.ews_to_abundances).
#
# A quick disclaimer: this example is intended to demonstrate the usage of Korg's fitting
# functionality, not to be an ironclad spectral analysis.
#
# We'll use a few packages in addition to Korg in this example. If you don't have them installed, you
# can run `using Pkg; Pkg.add(["CSV", "DataFrames", "PyPlot"])` to install them.

using Suppressor # hide #md
@suppress begin # hide #md
    using Korg, PythonPlot, CSV, DataFrames
end # hide #md

# # Reading in the data
# First, we need to read in the linelist, "window list" (the locations of lines to be fit), and
# spectrum. You can download the data files from Korg's GitHub repository at
# [https://github.com/ajwheeler/Korg.jl/tree/v1.0.0/docs/src/assets/Griffith_2022](https://github.com/ajwheeler/Korg.jl/tree/v1.0.0/docs/src/assets/Griffith_2022).
#
# ## Linelist
# As in [the Equivalent Widths Example](@ref EW_linelist), we want to reproduce the results of a
# specific paper, so we'll use the linelist from that paper.

linetable = CSV.read("../../assets/Griffith_2022/lines.csv", DataFrame);
linelist = Korg.Line.(Korg.air_to_vacuum.(linetable.wave_A),
                      linetable.loggf,
                      Korg.Species.(linetable.element),
                      linetable.lower_state_eV,
                      linetable.rad,
                      linetable.stark,
                      linetable.waals)

# Next, we'll read in the "window list", which provides wavelength ranges for features to be fit.
# We'll read this into a dictionary that maps atomic number to a vector of (lower, upper) wavelength
# bounds.

windowtable = CSV.File("../../assets/Griffith_2022/windows.tsv"; delim='\t');
windows = Dict()
for row in windowtable
    #atomic number
    Z = Korg.get_atoms(Korg.Species(row.species))[1]
    #get the windows for this element so far
    wins = get(windows, Z, [])
    #add the new window
    push!(wins, (Korg.air_to_vacuum(row.wave_base * 10), Korg.air_to_vacuum(row.wave_top * 10)))
    #update the dictionary
    windows[Z] = wins
end

# Finally, read in the observed spectrum.

spec = CSV.read("../../assets/Griffith_2022/2MASS_J03443498+0553014.csv", DataFrame)
spec.waveobs = Korg.air_to_vacuum.(spec.waveobs * 10)

# # Fitting iron lines to get stellar parameters
#
# First, we'll provide an initial guess for each parameter we want to fit.

# this is a "NamedTuple", but you can also use a dictionary if you prefer
initial_guess = (; Teff=5400, logg=3.8, M_H=-1.1, vmic=1.0)

# Next, we'll fit the spectrum using only the iron line locations from `windows`.
# The object produced, `fit_result`, contains several bits of information, most importantly the best-fit
# parameters.

winds = windows[26] # use Fe windows
fit_result = Korg.Fit.fit_spectrum(spec.waveobs, spec.flux, spec.err, linelist, initial_guess;
                                   windows=winds, R=50_000)
fit_result.best_fit_params

# It also contains the trace, which we can plot to see how the fit converged.
# Here's the $\chi^2$ as a function of the optimizer step.

figure(; figsize=(3, 3)) # hide #md
plot([t["chi2"] for t in fit_result.trace])
ylabel(L"χ^2")
xlabel("optimizer step")
gcf() # hide #md

# Here's the temperature as a function of the optimizer step.

figure(; figsize=(3, 3)) # hide #md
plot([t["Teff"] for t in fit_result.trace])
ylabel(L"$T_\mathrm{eff}$ [K]")
xlabel("optimizer step")
gcf() # hide #md

# The `fit_result` object also contains the best-fit spectrum, which we can plot in comparison to
# the observed one. Let's look at one of the iron lines.

# get the observed spectrum only at the wavelengths within a fitting window
obs_wls = spec.waveobs[fit_result.obs_wl_mask]
obs_flux = spec.flux[fit_result.obs_wl_mask]
obs_err = spec.err[fit_result.obs_wl_mask]

w = rand(winds) # choose a window around a random Fe line

# create a bitmask to plot the window plus 1 Å on each side for context
mask = w[1] - 1 .< obs_wls .< w[2] + 1

figure(; figsize=(3, 3)) # hide #md
scatter(obs_wls[mask], fit_result.best_fit_flux[mask]; c="r", label="Korg")
errorbar(obs_wls[mask], obs_flux[mask]; yerr=obs_err[mask], ls="", c="k", label="data")
legend()
gcf() # hide #md

# # Fit individual abundances
#
# Now, let's fit for the sodium abundance, holding the stellar parameters fixed.
# We'll use the stellar parameters from Griffith, rather than the ones from the analysis above.
# Comparing the abundances you get using each is left as an exercise to the reader.
#
# First, we'll define the initial guess for the Na abundance and the parameters to hold fixed.

init_params = (; Na=-1.0)
griffith_params = (Teff=5456, logg=3.86, M_H=-1.22, vsini=2.4, vmic=1.23)

# Next, we'll fit the sodium lines.
# The only difference from how we called `fit_spectrum` the first time is that we pass in a
# second set of parameters to hold fixed. [`Korg.Fit.fit_spectrum`](@ref) supports any combination
# of parameters to fit and hold fixed.

winds = windows[11] # windows for Na (atomic number 11)
na_result = Korg.Fit.fit_spectrum(spec.waveobs, spec.flux, spec.err, linelist,
                                  init_params, griffith_params;
                                  R=50_000, windows=winds) # 11 is the atomic number of Na
na_result.best_fit_params["Na"]

# # Parameter uncertainty
#
# Under the hood, [`Korg.Fit.fit_spectrum`](@ref) uses
# [LBFGS](https://en.wikipedia.org/wiki/Limited-memory_BFGS) to find the best-fit parameters.
# That means it produces an estimate of the Hessian of the likelihood function, and thus the
# covariance matrix of the best-fit parameters.

fit_result.covariance

# We caution the user that this is a very rough estimate, likely
# appropriate only for identifying pathological cases or the order of magnitude of the uncertainties.
