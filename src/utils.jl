using Statistics: quantile

normal_pdf(Δ, σ) = exp(-0.5*Δ^2 / σ^2) / √(2π) / σ

"""
    constant_R_LSF(flux, wls, R)

Applies a gaussian line spread function the the spectrum with flux vector `flux` and wavelength
vector `wls` with constant spectral resolution, R = λ/Δλ.

This will have weird behavior if your wavelength grid is not locally linearly-spaced.
It is intended to be run on a fine wavelength grid (``\\Delta\\lambda \\lesssim 0.05 \\AA``), then 
downsampled to the observational (or otherwise desired) grid.
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
Rectify the spectrum with flux vector `flux` and wavelengths `wls` by dividing out a moving 
`q`-quantile with `bandwidth`.
"""
function rectify(flux::AbstractVector{F}, wls; bandwidth=50, q=0.95) where F <: Real
    lb = 1
    ub = 1
    moving_mean = map(wls) do λ
        #move_bounds is defined in line_opacity.jl
        lb, ub = move_bounds(wls, lb, ub, λ, bandwidth)
        quantile(flux[lb:ub], q)
    end
    flux ./ moving_mean
end

"""
    air_to_vacuum(λ)

convert λ from an air wavelength to a vacuum wavelength
"""
function air_to_vacuum(λ)
    s = 1e4/λ
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s^2) + 0.0001599740894897 / (38.92568793293 - s^2)
    λ * n
end

"""
    vacuum_to_air(λ)

convert λ from a vacuum wavelength to an air wavelength
"""
function vacuum_to_air(λ)
    s = 1e4/λ
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s^2) + 0.00015998 / (38.9 - s^2)
    λ / n
end
