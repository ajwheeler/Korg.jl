normal_pdf(Δ, σ) = exp(-0.5*Δ^2 / σ^2) / √(2π) / σ

"""
Applies a gaussian line spread function with constant spectral resolution, R = λ/Δλ.

This will have weird behavior if your wavelength grid is not locally linearly-spaced.
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

function air_to_vacuum(λ)
    s = 1e4/λ
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s^2) + 0.0001599740894897 / (38.92568793293 - s^2)
    λ * n
end

function vacuum_to_air(λ)
    s = 1e4/λ
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s^2) + 0.00015998 / (38.9 - s^2)
    λ / n
end
