# [Tutorial notebooks](@id tutorials)
Here are a few tutorials in the form of Jupyter Notebooks to help you get started using Korg.  Each of the following bullets link to a folder on GitHub containing the associated Julia and Python notebooks (the `.ipynb` files) and any supporting data.  You can download everything if you want to run the notebook yourself, or just view the notebook on GitHub.
- [The Basics](https://github.com/ajwheeler/Korg.jl/tree/main/misc/Tutorial%20notebooks/basics) demonstrates how to read linelists, read or interpolate model atmospheres, and synthesize spectra. It also goes into some of the other data returned when you synthesize a spectrum.
- [Fitting](https://github.com/ajwheeler/Korg.jl/tree/main/misc/Tutorial%20notebooks/fitting) demonstrates the use of [`Korg.Fit.fit_spectrum`](@ref) to fit observational data.
- [Fitting with equivalent widths](https://github.com/ajwheeler/Korg.jl/tree/main/misc/Tutorial%20notebooks/EW%20fitting) demonstrates how to fit equivalent widths.

# Slurm example
[There is an example script](https://github.com/ajwheeler/Korg.jl/tree/main/misc/examples/generate_grid_with_slurm.jl) that shows how you might use Korg to generate many spectra within a Slurm job. Note that this example creates a file for each spectrum, which may not be desireable on your HPC environment.

# APOGEE example
[Here is an exmple script](https://github.com/ajwheeler/Korg.jl/tree/main/misc/examples/synthesize_apogee_spectrum.jl) showing how you might synthsize an APOGEE spectrum.  It demonstrates some bells and whistles: a precomputed LSF matrix and precomputed molecular cross sections.
