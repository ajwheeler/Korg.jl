using Statistics: quantile
using Interpolations: LinearInterpolation

normal_pdf(Δ, σ) = exp(-0.5*Δ^2 / σ^2) / √(2π) / σ

"""
    constant_R_LSF(flux, wls, R)

Applies a gaussian line spread function the the spectrum with flux vector `flux` and wavelength
vector `wls` with constant spectral resolution, ``R = \\lambda/\\Delta\\lambda.``

!!! warning
    This is a naive, slow implementation.  Do not use it when performance matters.

    `constant_R_LSF` will have weird behavior if your wavelength grid is not locally linearly-spaced.
    It is intended to be run on a fine wavelength grid, then downsampled to the observational (or 
    otherwise desired) grid.
"""
function constant_R_LSF(flux::AbstractVector{F}, wls, R) where F <: Real
    #ideas - require wls to be a range object? Use erf to account for grid edges?
    convF = zeros(F, length(flux))
    for i in 1:length(wls)
        λ0 = wls[i]
        σ = λ0 / R / 2
        mask = λ0 - 4σ .< wls .< λ0 + 4σ
        ϕ = normal_pdf.(wls[mask] .- λ0, λ0 / R / 2)
        ϕ ./= sum(ϕ)
        convF[mask] += flux[i]*ϕ
    end
    convF
end

"""
    rectify(flux, wls; bandwidth=50, q=0.95, wl_step=1.0)

Rectify the spectrum with flux vector `flux` and wavelengths `wls` by dividing out a moving
`q`-quantile with window size `bandwidth`.  `wl_step` controls the size of the grid that the moving 
quantile is calculated on and interpolated from.  Setting `wl_step` to 0 results in the exact 
calculation with no interpolation, but note that this is very slow.  

Experiments on real spectra show an agreement between the interpolated rectified spectrum and the 
"exact" one (with default values) at the 3 × 10^-4 level.
"""
function rectify(flux::AbstractVector{F}, wls; bandwidth=50, q=0.95, wl_step=1.0) where F <: Real
    inds = 1 : max(1, Int(floor(wl_step/step(wls)))) : length(wls)
    lb = 1
    ub = 1
    moving_quantile = map(wls[inds]) do λ
        #move_bounds is defined in line_opacity.jl
        lb, ub = move_bounds(wls, lb, ub, λ, bandwidth)
        quantile(flux[lb:ub], q)
    end
    itp = LinearInterpolation(wls[inds], moving_quantile)
    flux ./ itp.(wls)
end

"""
    air_to_vacuum(λ; cgs=λ<1)

Convert λ from an air wavelength to a vacuum wavelength.  
λ is assumed to be in Å if it is ⩾ 1, in cm otherwise.
Formula from Birch and Downs (1994) via the VALD website.
"""
function air_to_vacuum(λ; cgs=λ<1)
    if cgs
        λ *= 1e8
    end
    s = 1e4/λ
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s^2) + 0.0001599740894897 / (38.92568793293 - s^2)
    if cgs
        λ *= 1e-8
    end
    λ * n
end

"""
    vacuum_to_air(λ; cgs=λ<1)

convert λ from a vacuum wavelength to an air wavelength
λ is assumed to be in Å if it is ⩾ 1, in cm otherwise.
Formula from Birch and Downs (1994) via the VALD website.
"""
function vacuum_to_air(λ; cgs=λ<1)
    if cgs
        λ *= 1e8
    end
    s = 1e4/λ
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s^2) + 0.00015998 / (38.9 - s^2)
    if cgs
        λ *= 1e-8
    end
    λ / n
end

"""
    parse_fwf(source, rowspec; datarow, lastrow)

Parse fixed-width files, returning a vector of `NamedTuples`.
`source` can be either an `AbstractVector` of `Strings` or a filename.
`rowspec` should be a vector of tuples each containing, in order
- a range object corresponding to the column postion
- the type.  Strings will be pulled directly & other types will be parsed with `parse`.
- the column name, a `Symbol`.
- optionally, a function to apply to the parsed value.
"""
function parse_fwf(filename::String, rowspec; datarow=1, lastrow=0)
    open(filename, "r") do file
        for _ in 1:(datarow-1)
            readline(file)
        end
        _parse_fwf_core(eachline(file), rowspec, lastrow-datarow+1)
    end
end
function parse_fwf(lines::AbstractVector{String}, rowspec; datarow=1, lastrow=0)
    lastrow = lastrow == 0 ? length(lines) : lastrow
    _parse_fwf_core(lines[datarow:lastrow], rowspec)
end

function _parse_fwf_core(lines, rowspec, n=0)
    #use n
    map(lines) do line
        #broadcasts over each itemspec and passes to NamedTuple constuctor as kwargs
        (; _parse_item.(line, rowspec)...)
    end
end

function _parse_item(line, itemspec)
    if length(itemspec) == 4
        r, t, n, f = itemspec
    else
        r, t, n = itemspec
        f = identity
    end
    if t == String
        (n, f(line[r]))
    elseif r[2] > length(line) || strip(line[r]) == ""
        (n, zero(t))
    else
        (n, f(parse(t, line[r])))
    end
end
