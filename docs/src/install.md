## [1. Install Julia](@id install)
Install Julia by downloading a binary from [the website](https://julialang.org/downloads/).  
See also [the platform specific instructions](https://julialang.org/downloads/platform/), especially
if you would like to use Julia from the command line.  

!!! warning
    The latest versions of Korg support Julia 1.7 and greater, so make sure you are not using an 
    older version.  (Versions can be installed alongside eachother, and the 1.X versions are 
    backwards compatible.)

## 2. Install Korg

1. Lauch a julia session (either by launching the app, or typing `julia` on the command line if you have that set up).
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
    notebooks), and [PyPlot](https://github.com/JuliaPy/PyPlot.jl) (for calling `matplotlib` from 
    Julia). 

    If you want PyPlot to use you existing python/matplotlib installation, just do 
    `ENV["PYTHON"] = "/path/to/python"` before you install.  See [here](https://github.com/JuliaPy/PyCall.jl#specifying-the-python-version)
    for details.

## 3. (Optional) setup `PyJulia`.
We recommend using Korg from Julia, but if you would like to use Korg from Python, you can use 
[`PyJulia`](https://pyjulia.readthedocs.io/en/latest/index.html).  Their documentation has [detailed 
installation notes](https://pyjulia.readthedocs.io/en/latest/installation.html), but here's the 
short version: 
```
$ python3 -m pip install --user julia
$ python
>>> import julia
>>> julia.install()
```
