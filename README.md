[![Documentation](https://img.shields.io/badge/Documentation-blue.svg)](https://ajwheeler.github.io/Korg.jl/stable/)
[![Tests](https://github.com/ajwheeler/Korg.jl/actions/workflows/Test.yml/badge.svg)](https://github.com/ajwheeler/Korg.jl/actions/workflows/Test.yml)
[![codecov](https://codecov.io/gh/ajwheeler/Korg.jl/branch/main/graph/badge.svg?token=XXK2G8T8CJ)](https://codecov.io/gh/ajwheeler/Korg.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Code papers (please cite these if you use Korg): 
- [Korg: A Modern 1D LTE Spectral Synthesis Package](https://ui.adsabs.harvard.edu/abs/2023AJ....165...68W/abstract)
- [Korg: fitting, model atmosphere interpolation, and Brackett lines](https://ui.adsabs.harvard.edu/abs/2023arXiv231019823W/abstract) 

Tutorials:
- [Korg basics](https://github.com/ajwheeler/Korg.jl/blob/main/misc/Tutorial%20notebooks/basics/Basics.ipynb)
- [Korg basics in python](https://github.com/ajwheeler/Korg.jl/blob/main/misc/Tutorial%20notebooks/basics/Python%20Basics.ipynb)
- [other tutorial notebooks](https://github.com/ajwheeler/Korg.jl/tree/main/misc/Tutorial%20notebooks)


## Example
```julia
using Korg, PyPlot

wls, flux, continuum = synth(
    Teff=5000,
    logg=4.32,
    m_H=-1.1,
    C=-0.5,
    linelist=Korg.get_GALAH_DR3_linelist(),
    wavelengths=(5850, 5900)
)

# plot
figure(figsize=(12, 4))
plot(wls, flux, "k-")
xlabel(L"$\lambda$ [Å]")
ylabel(L"$F_\lambda/R_\mathrm{star}^2$ [erg s$^{-1}$ cm$^{-5}$]");
```

![image](https://github.com/ajwheeler/Korg.jl/assets/711963/70a13b45-4db2-472c-9121-fdd818a47105)


## You can also call Korg from python
See [the documentation](https://ajwheeler.github.io/Korg.jl/stable/install/) for setup instructions.
```python
from juliacall import Main as jl
jl.seval("using Korg"); Korg = jl.Korg

wls, flux, continuum = Korg.synth(
    Teff=5000,
    logg=4.32,
    m_H=-1.1,
    C=-0.5,
    linelist=Korg.get_GALAH_DR3_linelist(),
    wavelengths=(5850, 5900)
)
```
