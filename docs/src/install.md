While Korg is written in Julia, it's possible to use it from either Julia or Python.

# Using Korg from Julia

## [1. Install Julia](@id install)
Install Julia, by following the [instructions on the website](https://julialang.org/downloads/).

!!! warning
    The latest versions of Korg support Julia 1.7 and greater, so make sure you are not using an
    older version.  (Versions can be installed alongside eachother, and the 1.X versions are
    backwards compatible.)

## 2. Install Korg

1. Launch a julia session (either by launching the app, or typing `julia` on the command line if you have that set up).
2. Type `]` to enter `Pkg` mode.
3. Type `add Korg` to install Korg and its dependencies.
4. Press backspace or `CTRL+C` to exit `Pkg` mode (and return to the Julia REPL)

Alternatively, you can run
```julia
julia> using Pkg
julia> Pkg.add("Korg")
```

!!! tip
    If you are coming from Python, we also recommend installing
    [IJulia](https://github.com/JuliaLang/IJulia.jl) (for using Julia from Jupyter/IPython
    notebooks), and [PythonPlot](https://github.com/JuliaPy/PythonPlot.jl) (for calling `matplotlib` from Julia).

# Using Korg from Python
The recommended way to call Korg from Python is to use [juliacall](https://pypi.org/project/juliacall/).
Here's the quick version:

Install juliacall with
```bash
pip install juliacall
pip install juliapkg
```
Then, to install Korg, do this from a Python shell:
```python
import juliapkg
juliapkg.add("Korg", "acafc109-a718-429c-b0e5-afd7f8c7ae46")
juliapkg.resolve()
```
That's it! To use Korg from Python, just put these lines at the top of your script/notebook.
```python
from juliacall import Main as jl
jl.seval("using Korg")
Korg = jl.Korg
```

# Keeping Korg updated
In order to update Korg in the future, you can type:

Julia:
```julia
using Pkg
Pkg.update("Korg")
```

Python:
```python
from juliacall import Main as jl
jl.seval("using Pkg")
jl.Pkg.update("Korg")
```
