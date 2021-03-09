"""
    line_profile(temp, atomic_mass, line, wl)

The line profile, ϕ, at wavelengths `wls` in Ångstroms.
`temp` should be K, `atomic_mass` should be in g.
`line` should be one of the entries returned by `read_line_list`.
"""
function line_profile(temp::F, atomic_mass::F, line::NamedTuple, wls::AbstractVector{F}
                     ) where F <:  AbstractFloat
    #work in cgs
    λs = wls .* 1e-8 
    λ0 = line.wl * 1e-8

    #doppler-broadening parameter
    Δλ_D = λ0 * sqrt(2kboltz_cgs*temp / atomic_mass) / c_cgs

    #get all damping params from line list.  There may be better sources for this.
    #TODO these are log 10, right?
    γ = 10^line.log_gamma_rad + 10^line.log_gamma_stark + 10^line.log_gamma_vdW

    #Voigt function parameters
    v = @. abs(λs-λ0) / Δλ_D
    a = @. λs^2/(4π * c_cgs) * γ / Δλ_D

    @. voigt(a, v) / sqrt(π) / Δλ_D * 1e-8 #convert back to Å^-1
end

function harris_series(v) # assume v < 5
    H₀ = exp(-(v^2))
    H₁ = if (v < 1.3)
        -1.12470432 - 0.15516677v + 3.28867591v^2 - 2.34357915v^3 + 0.42139162v^4
    elseif v < 2.4
        -4.48480194 + 9.39456063v - 6.61487486v^2 + 1.98919585v^3 - 0.22041650v^4
    else #v < 5
        (0.554153432 + 0.278711796v - 0.188325687v^2 + 0.042991293v^3 - 0.003278278v^4) / (v^2 - 3/2)
    end
    H₂ = (1 - 2v^2) * H₀
    H₀, H₁, H₂
end

voigt(x, σ, γ) = lbvoigt(γ/2σ, abs(x)/σ)
function voigt(α, v)
    if α <= 0.2 
        if (v >= 5)
            (α / sqrt(π) / v^2) * (1 + 3/(2v^2) + 15/(4v^4))
        else
            H₀, H₁, H₂ = harris_series(v)
            H₀ + H₁*α + H₂*α^2
        end
    else  #unusual: α > 0.2
        if (α <= 1.4) && (α + v < 3.2)
            #modified harris series: M_i is H'_i in the source text
            H₀, H₁, H₂ = harris_series(v)
            M₀ = H₀
            M₁ = H₁ + 2/sqrt(π) * M₀
            M₂ = H₂ - M₀ + 2/sqrt(π) * M₁
            M₃ = 2/(3sqrt(π))*(1 - H₂) - 2/3 * v^2 * M₁ + 2/sqrt(π) * M₂
            M₄ = 2/3 * v^4 * M₀ - 2/(3sqrt(π)) * M₁ + 2/sqrt(π) * M₃
            ψ = 0.979895023 - 0.962846325α + 0.532770573α^2 - 0.122727278α^3
            ψ * (M₀ + M₁*α + M₂*α^2 + M₃*α^3 + M₄*α^4)
        else #α > 1.4 or (α > 0.2 and α + v > 3.2)
            u = sqrt(2) * (α^2 + v^2)
            sqrt(2/π) * (α/u) * (1 + (3v^2 - α^2)/u^2 + (15v^4 - 30α^2*v^2 + 2α^4)/u^4)
        end
    end
end

