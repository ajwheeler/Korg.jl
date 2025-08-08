# [Download this page as a Jupyter notebook](./EWs.ipynb) #md
# You may have to execute these cells to see the plot outputs #nb

# # Fitting "via synthesis" (direct fitting)
#
# In this example, we'll show how to use `Korg.Fit.fit_spectrum` to fit one of the spectra from
# [Griffith et al. 2022](https://ui.adsabs.harvard.edu/abs/2022arXiv221001821G). You will probably
# find it helpful to look at the [documentation for this
# function](https://ajwheeler.github.io/Korg.jl/stable/API/#Korg.Fit.fit_spectrum) as well. For
# fitting equivalent widths (rather than spectra directly) see [the documentation for
# `Korg.Fit.ews_to_abundances`](https://ajwheeler.github.io/Korg.jl/stable/API/#Korg.Fit.ews_to_abundances).
#
# This example is intended to demonstrate the usages of Korg's fitting functionality, not an
# ironclad spectral analysis.
#
# We'll use a few packages in addition to Korg in this example.  If don't have them installed, you
# can run `using Pkg; Pkg.add(["CSV", "DataFrames", "PyPlot"])` to install them.

using Suppressor # hide #md
@suppress begin # hide #md
    using Korg, PythonPlot, CSV, DataFrames
end # hide #md

# ## Reading in the data
# First, we need to read in the linelist, "window list" (the locations of lines to be fit), and
# spectrum.
#
# For this example, we'll use the same linelist and windows as Griffith et al. 2022, which means
# reading the CSV files they are stored in and getting the data into the format expected by Korg.
# For other linelists options, see the documentation for [Korg's built-in
# linelists](https://ajwheeler.github.io/Korg.jl/stable/API/#Built-in-linelists), and
# [`Korg.read_linelist`](https://ajwheeler.github.io/Korg.jl/stable/API/#Korg.read_linelist).

# TODO
