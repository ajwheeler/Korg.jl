
## 1. Install Julia
Install Julia by downloading a binary from [the website](https://julialang.org/downloads/).  
See also [the platform specific instructions](https://julialang.org/downloads/platform/), especially
if you would like to use Julia from the command line.  

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
    Julia.


## 3. (Optional) setup `PyJulia`.
If you would like to use Korg from Python, you can use 
[`PyJulia`](https://pyjulia.readthedocs.io/en/latest/index.html).  Their documentation has [detailed 
installation notes](https://pyjulia.readthedocs.io/en/latest/installation.html), but here's the 
short version: 
```
$ python3 -m pip install --user julia
$ python
>>> import julia
>>> julia.install()
```
