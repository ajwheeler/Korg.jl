# [Download this page as a Jupyter notebook](./Equivalent Widths.ipynb) #md
# You may have to execute these cells to see the plot outputs #nb

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
    using Korg, PythonPlot, CSV, DataFrames
end # hide #md

# ## [Creating a custom linelist](@id EW_linelist)
#
# Korg can read in linelists in a variety of formats, including the VALD, MOOG, Kurucz, and ExoMol
# formats, but sometimes you'll want to use one in another format, e.g. from a table in a paper.
# We'll use the data from Table 1 of
# [Melendez 2014](https://ui.adsabs.harvard.edu/abs/2014ApJ...791...14M/abstract),
# which are available in the
# [`Table 1.dat` file](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJ/791/14#/browse).
#
# First, we read the data into a DataFrame called `lines`. This is a convenient way to work with
# tabular data.

lines = CSV.read("../../assets/Table 1.dat", DataFrame; skipto=25, delim=' ', ignorerepeated=true,
                 header=["wl", "species", "ExPot", "log_gf", "C6", "EW_18Sco", "EW_Sun"]);

# Next, we convert the numbers in the species column to [`Korg.Species`](@ref) objects.
# The `Korg.Species` constructor takes a string and returns a `Korg.Species` object, and supports
# nearly all of the formats for specifying species found in the wild.
#
# The `.` applies the function to each element of the series. This is called "broadcasting".

lines.species = Korg.Species.(lines.species)

# Let's look at Fe lines only, and sort them by wavelength.

filter!(lines) do row
    #careful, get_atom throws an error when applied to a molecular species
    Korg.get_atom(row.species) == 26 # atomic number of Fe
end

sort!(lines, :wl)
lines[1:4, :] # look at the first few rows

# To pass the lines to Korg, we need to turn each row of the `lines` DataFrame into a
# [`Korg.Line`](@ref) object. Korg will use reasonable defaults for the broadening parameters,
# but see the [`Korg.Line` documentation](@ref Korg.Line) for details on how to specify them if you need
# to do that.

linelist = Korg.Line.(lines.wl, # can be in either cm or Å (like these), but NOT nm
                      lines.log_gf,
                      lines.species, # needs to be a Korg.Species, which we handled in the cell above
                      lines.ExPot) # excitation potential, i.e. lower level energy (must be in eV)

# We can use `lines` directly with Korg at this point, but if we want to save it for later use,
# we can save it to an hdf5 file.

#read this back in with Korg.read_linelist("Fe_lines.h5"; format="korg")
Korg.save_linelist("Fe_lines.h5", linelist);

# ## Equivalent width fitting

# Now that we have our linelist, let's perform equivalent width fitting to determine abundances.
# We'll compare the Sun and 18 Sco using the equivalent widths from the paper.
#
# First, we'll define the parameters for the Sun and 18 Sco.

#solar params
sun_Teff, sun_logg, sun_Fe_H, sun_vmic = 5777, 4.44, 0.0, 1.0
#vector of abundances for the sun
sun_A_X = Korg.format_A_X(sun_Fe_H)
#model atmosphere for the sun
sun_atm = Korg.interpolate_marcs(sun_Teff, sun_logg, sun_A_X)
# and likewise for 18 Sco
sco_teff, sco_logg, sco_fe_h, sco_vmic = (5823, 4.45, 0.054, sun_vmic + 0.02)
sco_A_X = Korg.format_A_X(sco_fe_h)
sco_atm = Korg.interpolate_marcs(sco_teff, sco_logg, sco_A_X)

# Now we can calculate abundances from the EWs for each star.

A_sun = Korg.Fit.ews_to_abundances(sun_atm, linelist, sun_A_X, lines.EW_Sun; vmic=sun_vmic)
A_18Sco = Korg.Fit.ews_to_abundances(sco_atm, linelist, sco_A_X, lines.EW_18Sco; vmic=sco_vmic);

# Let's plot the abundance differences as a function of excitation potential to check for
# non-LTE effects or other systematic issues.

neutrals = [spec.charge == 0 for spec in lines.species] #bitmask for the lines of Fe I vs Fe II
χ = [l.E_lower for l in linelist] # the excitation potential for each line

figure(; figsize=(8, 2)) # hide #md
scatter(χ[neutrals], (A_18Sco-A_sun)[neutrals]; label="Fe I")
scatter(χ[.!neutrals], (A_18Sco-A_sun)[.!neutrals]; label="Fe II")
ylabel(L"A(Fe)_\mathrm{18 Sco} - A(Fe)_\mathrm{sun}")
xlabel("χ [eV]")
legend()
gcf() # hide #md
