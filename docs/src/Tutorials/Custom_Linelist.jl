# [Download this page as a Jupyter notebook](./Custom_Linelist.ipynb) #md
# You may have to execute these cells to see the plot outputs TODO #nb

# # [Custom linelists](@id custom_linelists)
#
# Korg can read in linelists in a variety of formats, including the VALD, MOOG, Kurucz, and ExoMol
# formats, but sometimes you'll want to use one in another format, e.g. a custom CSV.
# This tutorial demonstrates how to read in custom linelist and save it in a format that Korg can
# easily use.
#
# First, we'll import the necessary packages.
#
# We are using the CSV and DataFrames packages to read in and process the linelist.

using Suppressor # hide #md
@suppress begin # hide #md
    using Korg, CSV, DataFrames
end # hide #md

# Let's use the linelist from [Melendez 2014](https://ui.adsabs.harvard.edu/abs/2014ApJ...791...14M/abstract).
# We'll use the data from Table 1,
# [which are available in the `Table 1.dat` file](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJ/791/14#/browse).
#
# First, we read the data into a DataFrame called "lines". This is a convienient way to work with
# tabular data.

lines = CSV.read("../../assets/Table 1.dat", DataFrame; skipto=25, delim=' ', ignorerepeated=true,
                 header=["wl", "species", "ExPot", "log_gf", "C6", "EW_18Sco", "EW_Sun"])

# Next, we convert the numbers in the species column to [`Korg.Species`](@ref) objects.
# The `Korg.Species` constructor takes a string, and returns a `Korg.Species` object, and supports
# nearly all of the formats for specifying species found in the wild.
#
# The `.` applies the function to each element of the series. This is called "broadcasting".

lines.species = Korg.Species.(lines.species)

# Let's look at Fe lines only, and sort them by wavelength.

filter!(lines) do row
    Korg.get_atoms(row.species) == [26]
end

sort!(lines, :wl)

# To pass the lines to Korg, we need to turn each row of the `lines` DataFrame into a
# [`Korg.Line`](@ref) object. This table doesn't contain any broadening parameters, so Korg will
# use reasonable defaults, but see the [`Korg.Line` documentation](@ref) details on how to specify
# that info if you have it.
linelist = Korg.Line.(lines.wl, # can be in either cm or Å (like these), but NOT nm
                      lines.log_gf,
                      lines.species, # needs to be a Korg.Species, which we handled in the cell above
                      lines.ExPot) # excitation potential, i.e. lower level energy (must be in eV)

# We can use `lines` directly with Korg at this point, but if we want to save it for later use,
# we can save it to an hdf5 file.

Korg.save_linelist("Fe_lines.h5", linelist)
;
