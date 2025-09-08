# Korg.jl

[![Documentation](https://img.shields.io/badge/Documentation-blue.svg)](https://ajwheeler.github.io/Korg.jl/stable/)
[![Tests](https://github.com/ajwheeler/Korg.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ajwheeler/Korg.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/ajwheeler/Korg.jl/branch/main/graph/badge.svg?token=XXK2G8T8CJ)](https://codecov.io/gh/ajwheeler/Korg.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Korg is a package for computing stellar spectra from 1D model atmospheres and linelists assuming local thermodynamic equilibrium. It can be used with Julia or Python. Here's some things it can do:
- Computing spectra from Teff, logg, abundances, etc.
- Fitting whole spectra or individual lines via synthesis or equivalent widths
- Excitation-ionization balance from equivalent width data
- Model atmosphere interpolation, parsing and using atmosphere files (MARCS)
- Parsing and using linelists in VALD, Kurucz, MOOG, ExoMol, and Turbospectrum formats, with several defaults built-in.
- Automatic differentiaion (via ForwardDiff.jl)
- Synthesis with arbitrary abundances/solar abundance scales, with several defaults built-in.

## Example
(Python version below)
```julia
using Korg, PythonPlot

wls, flux, continuum = synth(
    Teff=5000, # effective temperature of 5000 Kelvin
    logg=4.32, # surface gravity of 10^(4.32) cm/s²
    M_H=-1.1,  # metallicity, [M/H]. Overridden for individual elements by alpha_H and individual abundances
    C=-0.5,    # The Carbon abundance, [C/H].  Works for anything from He to U.
    linelist=Korg.get_GALAH_DR3_linelist(),
    wavelengths=(5850, 5900)
)

# plot
figure(figsize=(12, 4))
plot(wls, flux, "k-")
xlabel(L"$\lambda$ [Å]")
ylabel("continuum-normalized flux");
```
![spectrum](https://github.com/ajwheeler/Korg.jl/assets/711963/70a13b45-4db2-472c-9121-fdd818a47105)

See the [documentation for `synth`](https://ajwheeler.github.io/Korg.jl/stable/API/#Korg.synth), or [the documentation for `synthesize`](https://ajwheeler.github.io/Korg.jl/stable/API/#Korg.synthesize) for advanced usage.

# Code papers (please cite these if you use Korg):
- [Korg: A Modern 1D LTE Spectral Synthesis Package](https://ui.adsabs.harvard.edu/abs/2023AJ....165...68W/abstract). This is also a good overview of how spectral synthesis works, the inputs and outputs, etc.
- [Korg: fitting, model atmosphere interpolation, and Brackett lines](https://ui.adsabs.harvard.edu/abs/2023arXiv231019823W/abstract)

# Getting help
If you are having trouble using or installing Korg, please get in touch by [opening a GitHub issue](https://github.com/ajwheeler/Korg.jl/issues) (preferred), or [sending Adam an email](mailto:adamwhlr@gmail.com).

# You can also call Korg from python
See [the documentation](https://ajwheeler.github.io/Korg.jl/stable/install/) for setup instructions.
```python
from juliacall import Main as jl
jl.seval("using Korg"); Korg = jl.Korg

# calling Korg.synth is exactly the same as in Julia.
wls, flux, continuum = Korg.synth(
    Teff=5000,
    logg=4.32,
    M_H=-1.1,
    C=-0.5,
    linelist=Korg.get_GALAH_DR3_linelist(),
    wavelengths=(5850, 5900)
)
```

# Multithreading
Korg can use multithreading to speed up line opacity calculation, the most expensive step for syntheses.
To use it [launch Julia with more than one thread](https://docs.julialang.org/en/v1/manual/multi-threading/), using the `-t` command-line argument, or by setting the [`$JULIA_NUM_THREADS`](https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_NUM_THREADS) environment variable.
