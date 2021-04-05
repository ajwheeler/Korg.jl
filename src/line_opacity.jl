"""
    line_opacity(linelist, wls, temp, n_densities, atomic_masses, partition_fns, 
                 ionization_energies; window_size)

Calculate the opacity coefficient, α, in units of cm^-1 from all lines in `linelist`, at wavelengths
`wls`. 

other arguments:
- `temp` the temerature in K
- `n_densities`, a Dict mapping species to absolute number density [cm^-3].
- `window_size` (optional, default: 20), the maximum distance from the line center at which line 
opacities should be calculated in included.
- `ξ` is the microturbulent velocity in cm/s
"""
function line_absorption(linelist, wls, temp, n_densities::Dict, atomic_masses::Dict, 
                         partition_fns::Dict, ionization_energies::Dict, ξ
                         ; window_size=20.0)

    α_lines = zeros(length(wls))
    #lb and ub are the indices to the upper and lower wavelengths in the "window", i.e. the shortest
    #and longest wavelengths which feel the effect of each line 
    lb = 1
    ub = 1
    for line in linelist
        while wls[lb] < line.wl - window_size
            lb += 1
        end
        while wls[ub+1] < line.wl + window_size && ub+1 < length(wls)
            ub += 1
        end

        #line profile (normalized)
        ϕ = line_profile(temp, atomic_masses[get_elem(line.species)], ξ, line, view(wls,lb:ub))

        #cross section
        σ = sigma_line(line.wl, line.log_gf)

        #stat mech quantities
        #number density of particles in the relevant excitation state
        boltzmann_factor = exp(- line.E_lower / kboltz_eV / temp)
        n = n_densities[line.species] * boltzmann_factor / partition_fns[line.species](temp)
        #the factor (1 - exp(hν₀ / kT))
        levels_factor = 1 - exp(-c_cgs * hplanck_cgs / (line.wl * 1e-8) / kboltz_cgs / temp)

        @. α_lines[lb:ub] += ϕ * σ * n * levels_factor 
    end
    α_lines
end

"""
    sigma_line(wl, log_gf)

The cross-section at wavelength `wl` in Ångstroms of a transition for which the product of the
degeneracy and oscillator strength is `10^log_gf`.
"""
function sigma_line(wl, log_gf) where F <: AbstractFloat
    #work in cgs
    λ = wl .* 1e-8 
    e  = electron_charge_cgs
    mₑ = electron_mass_cgs
    c  = c_cgs

    (π*e^2/mₑ/c) * (λ^2/c) * 10^log_gf 
end

"""
    line_profile(temp, atomic_mass, ξ, line, wls)

The line profile, ϕ, at wavelengths `wls` in Ångstroms.
`temp` should be in K, `atomic_mass` should be in g, microturbulent velocity, `ξ`, should be in 
cm/s, `line` should be one of the entries returned by `read_line_list`.
Note that this returns values in units of cm^-1, not Å^-1
"""
function line_profile(temp::F, atomic_mass::F, ξ::F, line::NamedTuple, wls::AbstractVector{F}
                     ) where F <:  AbstractFloat
    #work in cgs
    λs = wls .* 1e-8 
    λ0 = line.wl * 1e-8

    #doppler-broadening parameter
    Δλ_D = λ0 * sqrt(2kboltz_cgs*temp / atomic_mass + ξ^2) / c_cgs

    #get all damping params from line list.  There may be better sources for this.
    γ = 10^line.log_gamma_rad + 10^line.log_gamma_stark + 10^line.log_gamma_vdW

    #Voigt function parameters
    v = @. abs(λs-λ0) / Δλ_D
    a = @. λs^2/(4π * c_cgs) * γ / Δλ_D

    @. voigt(a, v) / sqrt(π) / Δλ_D 
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
