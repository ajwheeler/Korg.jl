[![Tests](https://github.com/ajwheeler/Korg.jl/actions/workflows/Test.yml/badge.svg)](https://github.com/ajwheeler/Korg.jl/actions/workflows/Test.yml)
[![codecov](https://codecov.io/gh/ajwheeler/Korg.jl/branch/main/graph/badge.svg?token=XXK2G8T8CJ)](https://codecov.io/gh/ajwheeler/Korg.jl)

1D LTE stellar spectral synthesis for FGK stars in pure Julia.

[Code paper](https://ui.adsabs.harvard.edu/abs/2023AJ....165...68W/abstract) (please cite this if you use Korg)

[Tutorial](https://github.com/ajwheeler/Korg.jl/blob/main/misc/Tutorial%20notebooks/Tutorial.ipynb)
[(Python version)](https://github.com/ajwheeler/Korg.jl/blob/main/misc/Tutorial%20notebooks/Python%20Tutorial.ipynb)

[Documentation](https://ajwheeler.github.io/Korg.jl/stable/)

## Example
```julia
using Korg, PyPlot
lines = read_linelist("linelist.vald", format="vald")
atm = read_model_atmosphere("s6000_g+1.0_m0.5_t05_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
sol = synthesize(atm, lines, format_A_X(0), 5000, 5030);

figure(figsize=(12, 4))
plot(sol.wavelengths, sol.flux, "k-")
xlabel(L"$\lambda$ [Å]")
ylabel(L"$F_\lambda/R_\mathrm{star}^2$ [erg s$^{-1}$ cm$^{-5}$]");
```

![image](https://user-images.githubusercontent.com/711963/199083747-9d9d89b4-10a5-42f7-9273-11e9f6d2dfa1.png)

## You can also call Korg from python
See [the documentation](https://ajwheeler.github.io/Korg.jl/stable/install/) for setup instructions.
```python
import matplotlib.pyplot as plt
from julia import Korg

lines = Korg.read_linelist("linelist.vald", format="vald")
atm = Korg.read_model_atmosphere("s6000_g+1.0_m0.5_t05_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
A_X = Korg.format_A_X(0)
sol = Korg.synthesize(atm, lines, A_X, 5000, 5030)

fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(sol.wavelengths, sol.flux, "k-")
ax.set_xlabel("$\lambda$ [Å]")
ax.set_ylabel("$F_\lambda/R_\mathrm{star}^2$ [erg s$^{-1}$ cm$^{-5}$]")
```
