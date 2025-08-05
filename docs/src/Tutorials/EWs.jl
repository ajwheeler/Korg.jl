# [Download this page as a Jupyter notebook](./EWs.ipynb) #md
# You may have to execute these cells to see the plot outputs TODO #nb

# # Equivalent width fitting

# This tutorial shows how to use Korg to fit the abundances of a star using equivalent widths.
# We'll use the data from [Melendez 2014](https://ui.adsabs.harvard.edu/abs/2014ApJ...791...14M/abstract),
# which are available in the `Table 1.dat` file, and recreate some of the analysis from that paper.
#
# First, we'll import the necessary packages.
# These examples use the PythonPlot package, which provides a nice Julia interface to Matplotlib,
# but of course you can plot however you like.
# We'll also use CSV and DataFrames to read in the data from the paper.

using Suppressor # hide #md
@suppress begin # hide #md
    using Korg, PyPlot, CSV, DataFrames
end # hide #md

# First, read in the data from the paper.

lines = CSV.read("../../assets/Table 1.dat", DataFrame; skipto=25, delim=' ', ignorerepeated=true,
                 header=["wl", "species", "ExPot", "log_gf", "C6", "EW_18Sco", "EW_Sun"])

# The [Custom linelists](@ref custom_linelists) example shows how to turn the line data in this
# table into a linelist that Korg can use. Let's skip that step for now, and just read in the
# linelist from the hdf5 file we wrote in the [Custom linelists](@ref custom_linelists) example.

linelist = Korg.read_linelist("Fe_lines.h5"; format="korg")

# TODO

# solar params
sun_Teff, sun_logg, sun_Fe_H, sun_vmic = 5777, 4.44, 0.0, 1.0

# vector of abundances for the sun
sun_A_X = Korg.format_A_X(sun_Fe_H)

# interpolate a model atmosphere for the sun
sun_atm = Korg.interpolate_marcs(sun_Teff, sun_logg, sun_A_X)

# and likewise for 18 Sco
sco_teff, sco_logg, sco_fe_h, sco_vmic = (5823, 4.45, 0.054, sun_vmic + 0.02)
sco_A_X = Korg.format_A_X(sco_fe_h)
sco_atm = Korg.interpolate_marcs(sco_teff, sco_logg, sco_A_X)

# calculate abundances from the EWs for each star
# TODO

@time A_sun = Korg.Fit.ews_to_abundances(sun_atm, linelist, sun_A_X, lines.EW_Sun;
                                         vmic=sun_vmic)
@time A_18Sco = Korg.Fit.ews_to_abundances(sco_atm, linelist, sco_A_X, lines.EW_18Sco;
                                           vmic=sco_vmic)

# TODO
# get a bitmask for the lines of Fe I vs Fe II
neutrals = [spec.charge == 0 for spec in lines.species]
χ = [l.E_lower for l in linelist] # the excitation potential for each line

figure(; figsize=(8, 2)) # hide #md
scatter(χ[neutrals], (A_18Sco-A_sun)[neutrals]; label="Fe I")
scatter(χ[.!neutrals], (A_18Sco-A_sun)[.!neutrals]; label="Fe II")
ylabel(L"A(Fe)_\mathrm{18 Sco} - A(Fe)_\mathrm{sun}")
xlabel("χ [eV]")
legend()
gcf() # hide #md
;
