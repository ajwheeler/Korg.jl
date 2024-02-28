[![Tests](https://github.com/ajwheeler/Korg.jl/actions/workflows/Test.yml/badge.svg)](https://github.com/ajwheeler/Korg.jl/actions/workflows/Test.yml)
[![codecov](https://codecov.io/gh/ajwheeler/Korg.jl/branch/main/graph/badge.svg?token=XXK2G8T8CJ)](https://codecov.io/gh/ajwheeler/Korg.jl)

[Documentation](https://ajwheeler.github.io/Korg.jl/stable/)

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
lines = Korg.get_GALAH_DR3_linelist()
A_X = format_A_X(-1.1, Dict("C"=>-0.5))
atm = interpolate_marcs(5000.0, 4.32, A_X)
sol = synthesize(atm, lines, A_X, 5850, 5900);

figure(figsize=(12, 4))
plot(sol.wavelengths, sol.flux, "k-")
xlabel(L"$\lambda$ [Å]")
ylabel(L"$F_\lambda/R_\mathrm{star}^2$ [erg s$^{-1}$ cm$^{-5}$]");
```

![image](https://github.com/ajwheeler/Korg.jl/assets/711963/70a13b45-4db2-472c-9121-fdd818a47105)


## You can also call Korg from python
See [the documentation](https://ajwheeler.github.io/Korg.jl/stable/install/) for setup instructions.
```python
# imports
from juliacall import Main as jl
jl.seval("using Korg"); Korg = jl.Korg
from matplotlib import pyplot as plt

# synthesize spectrum
lines = Korg.get_GALAH_DR3_linelist()
A_X = Korg.format_A_X(-1.1, {"C": -0.5})
atm = Korg.interpolate_marcs(5000.0, 4.32, A_X)
sol = Korg.synthesize(atm, lines, A_X, 5850, 5900);

# plot
fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(sol.wavelengths, sol.flux, "k-")
ax.set_xlabel("$\lambda$ [Å]")
ax.set_ylabel("$F_\lambda/R_\mathrm{star}^2$ [erg s$^{-1}$ cm$^{-5}$]")
```
