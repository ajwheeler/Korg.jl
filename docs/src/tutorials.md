# [Other Tutorials](@id tutorials)

In addition to the examples here, there are a few more resources than may help you get started with
Korg.

## Slurm example
[There is an example script](https://github.com/ajwheeler/Korg.jl/tree/main/misc/examples/generate_grid_with_slurm.jl) that shows how you might use Korg to generate many spectra within a Slurm job. Note that this example creates a file for each spectrum, which may not be desireable on your HPC environment.

## APOGEE example
[Here is an exmple script](https://github.com/ajwheeler/Korg.jl/tree/main/misc/examples/synthesize_apogee_spectrum.jl) showing how you might synthsize an APOGEE spectrum.  It demonstrates some bells and whistles: a precomputed LSF matrix and precomputed molecular cross sections.
