[![Documentation](https://img.shields.io/badge/Documentation-blue.svg)](https://ajwheeler.github.io/Korg.jl/stable/)
[![Tests](https://github.com/ajwheeler/Korg.jl/actions/workflows/Test.yml/badge.svg)](https://github.com/ajwheeler/Korg.jl/actions/workflows/Test.yml)
[![codecov](https://codecov.io/gh/ajwheeler/Korg.jl/branch/main/graph/badge.svg?token=XXK2G8T8CJ)](https://codecov.io/gh/ajwheeler/Korg.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Code papers (please cite these if you use Korg):
- [Korg: A Modern 1D LTE Spectral Synthesis Package](https://ui.adsabs.harvard.edu/abs/2023AJ....165...68W/abstract). This is also a good overview of how spectral synthesis works, the inputs and outputs, etc.
- [Korg: fitting, model atmosphere interpolation, and Brackett lines](https://ui.adsabs.harvard.edu/abs/2023arXiv231019823W/abstract)

Tutorials:
- [Korg basics](https://github.com/ajwheeler/Korg.jl/blob/main/misc/Tutorial%20notebooks/basics/Basics.ipynb)
- [Korg basics in python](https://github.com/ajwheeler/Korg.jl/blob/main/misc/Tutorial%20notebooks/basics/Python%20Basics.ipynb)
- [other tutorial notebooks](https://github.com/ajwheeler/Korg.jl/tree/main/misc/Tutorial%20notebooks)


## Example
```julia
using Korg, PyPlot

wls, flux, continuum = synth(
    Teff=5000, # effective temperature of 5000 Kelvin
    logg=4.32, # surface gravity of 10^(4.32) cm/s²
    m_H=-1.1,  # metallicity, [m/H]. Overridden for individual elements by alpha_H and individual abundances
    C=-0.5,    # The Carbon abundance, [C/H].  Works for anything from He to U.
    linelist=Korg.get_GALAH_DR3_linelist(),
    wavelengths=(5850, 5900)
)

# plot
figure(figsize=(12, 4))
plot(wls, flux, "k-")
xlabel(L"$\lambda$ [Å]")
ylabel(L"$F_\lambda/R_\mathrm{star}^2$ [erg s$^{-1}$ cm$^{-5}$]");
```
See the [documentation for `synth`](https://ajwheeler.github.io/Korg.jl/stable/API/#Korg.synth), or [the documentation for `synthesize`](https://ajwheeler.github.io/Korg.jl/stable/API/#Korg.synth) for advanced usage.

![image](https://github.com/ajwheeler/Korg.jl/assets/711963/70a13b45-4db2-472c-9121-fdd818a47105)


## You can also call Korg from python
See [the documentation](https://ajwheeler.github.io/Korg.jl/stable/install/) for setup instructions.
```python
from juliacall import Main as jl
jl.seval("using Korg"); Korg = jl.Korg

# calling Korg.synth is exactly the same as in Julia.
wls, flux, continuum = Korg.synth(
    Teff=5000,
    logg=4.32,
    m_H=-1.1,
    C=-0.5,
    linelist=Korg.get_GALAH_DR3_linelist(),
    wavelengths=(5850, 5900)
)
```

# Multithreading
Korg can use multithreading to speed up line opacity calculation, the most epxensive step for for syntheses.
To use it [launch Julia with more than one thread](https://docs.julialang.org/en/v1/manual/multi-threading/), using the `-t` command-line argument, or by setting the [`$JULIA_NUM_THREADS`](https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_NUM_THREADS) environment variable.
