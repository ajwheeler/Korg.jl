normal_pdf(Δ, σ) = exp(-0.5*Δ^2 / σ^2) / √(2π) / σ

"""
Returns the matrix which transforms a spectrum to its LSF-convolved version for spectral resolution 
(λ / Δλ) R.

This will have weird behavior if your wavelength grid is not locally linearly-spaced.
It is intended to be run on a fine wavelength grid, then downsampled to the observational (or otherwise desired) grid.
"""
function line_spread_matrix(wls, R)
    M = Matrix{Float64}(undef, length(wls), length(wls))
    for i in 1:length(wls)
        λ0 = wls[i]
        M[i, :] = normal_pdf.(wls .- λ0, λ0 / R / 2)
        M[i, :] /= sum(M[i, :])
    end
    M
end
