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
