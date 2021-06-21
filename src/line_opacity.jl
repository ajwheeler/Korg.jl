using SpecialFunctions: gamma

"""
    line_absorption(linelist, λs, temp, nₑ, n_densities, partition_fns, ξ
                   ; α_cntm=nothing, cutoff_threshold=1e-3, window_size=20.0*1e-8)

Calculate the opacity coefficient, α, in units of cm^-1 from all lines in `linelist`, at wavelengths
`λs`. 

other arguments:
- `temp` the temerature in K
- `n_densities`, a Dict mapping species to absolute number density [cm^-3].
- `partition_fns`, a Dict containing the partition function of each species
- `ξ` is the microturbulent velocity in cm/s (n.b. NOT km/s)

The three keyword arguments specify how the size of the window over which each line is calculated 
is chosen.  
- If `α_cntm` is passed as a callable returning the continuum opacity as a function of 
wavelength, the line window will extend to 4 doppler widths or the wavelength at which the Lorentz
wings of the line are at `cutoff_threshold * α_cntm[line.wl]`, whichever is greater.  
`cutoff_threshold` defaults to 1e-3.  
- if `α_cntm` is not passed, defaults to `window_size`, which is 2e-7 (in cm, i.e. 20 Å) unless
otherwise specified
"""
function line_absorption(linelist, λs, temp, nₑ, n_densities::Dict, partition_fns::Dict, ξ 
                         ; α_cntm=nothing, cutoff_threshold=1e-3, window_size=20.0*1e-8)
    #type shenanigans to allow autodiff to do its thing
    α_type = typeof(promote(linelist[1].wl, λs[1], temp, nₑ, n_densities["H_I"], ξ, cutoff_threshold
                           )[1])
    α_lines = zeros(α_type, length(λs))

    if length(linelist) == 0
        return α_lines
    end

    #lb and ub are the indices to the upper and lower wavelengths in the "window", i.e. the shortest
    #and longest wavelengths which feel the effect of each line 
    lb = 1
    ub = 1
    for line in linelist
        mass = get_mass(strip_ionization(line.species))
        
        #doppler-broadening parameter
        Δλ_D = line.wl * sqrt(2kboltz_cgs*temp / mass + ξ^2) / c_cgs

        #get all damping params from line list.  There may be better sources for this.
        Γ = line.gamma_rad 
        if !ismolecule(line.species) 
            Γ += (nₑ*scaled_stark(line.gamma_stark, temp) + 
                  (n_densities["H_I"] + 0.42n_densities["He_I"])*scaled_vdW(line.vdW, mass, temp))
        end

        #doing this involves an implicit aproximation that λ(ν) is linear over the line window
        Δλ_L = Γ * line.wl^2 / c_cgs

        #cross section
        σ = sigma_line(line.wl)

        #stat mech quantities 
        #number density of particles in the relevant excitation state
        boltzmann_factor = exp(- line.E_lower / kboltz_eV / temp)
        #the factor (1 - exp(hν₀ / kT)) to account for stimulated emission
        emission_factor = 1 - exp(-c_cgs * hplanck_cgs / line.wl / kboltz_cgs / temp)
        levels_factor = boltzmann_factor / partition_fns[line.species](temp) * emission_factor

        line_amplitude = 10.0^line.log_gf * n_densities[line.species] * σ * levels_factor

        if !isnothing(α_cntm)
            α_crit = α_cntm(line.wl) * cutoff_threshold
            Δλ_crit = sqrt(line_amplitude * Δλ_L / α_crit) #where α from lorentz component == α_crit
            window_size = max(4Δλ_D, Δλ_crit)
        end
        lb, ub = move_bounds(λs, lb, ub, line.wl, window_size)
        if lb==ub
            continue
        end

        invΔλ_D = 1/Δλ_D
        @inbounds view(α_lines, lb:ub) .+= line_profile.(line.wl, invΔλ_D, Δλ_L, line_amplitude,
                                                         view(λs, lb:ub))
    end
    α_lines
end

"the stark broadening gamma scaled acording to its temperature dependence"
scaled_stark(γstark, T; T₀=10_000) = γstark * (T/T₀)^(1/6)

             
"""
the vdW broadening gamma scaled acording to its temperature dependence, using either simple scaling 
or ABO
"""
scaled_vdW(vdW::Real, m, T, T₀=10_000) = vdW * (T/T₀)^0.3
function scaled_vdW(vdW::Tuple{F, F}, m, T) where F <: Real
    v₀ = 1e6 #σ is given at 10_000 m/s = 10^6 cm/s
    σ = vdW[1]
    α = vdW[2]

    invμ = 1/(1.008*amu_cgs) + 1/m #inverse reduced mass
    vbar = sqrt(8 * kboltz_cgs * T / π * invμ) #relative velocity
    #n.b. "gamma" is the gamma function, not a broadening parameter
    2 * (4/π)^(α/2) * gamma((4-α)/2) * v₀ * σ * (vbar/v₀)^(1-α)
end


function move_bounds(λs::AbstractRange, lb, ub, λ₀, window_size)
    len = length(λs)
    lb = clamp(Int(cld(λ₀ - window_size - λs[1], step(λs)) + 1), 1, len)
    ub = clamp(Int(fld(λ₀ + window_size - λs[1], step(λs)) + 1), 1, len)
    lb,ub
end
#walk lb and ub to be window_size away from λ₀. assumes λs is sorted
function move_bounds(λs, lb, ub, λ₀, window_size)
    while lb+1 < length(λs) && λs[lb] < λ₀ - window_size
        lb += 1
    end
    while lb > 1 && λs[lb-1] > λ₀ - window_size
        lb -= 1
    end
    while ub < length(λs) && λs[ub+1] < λ₀ + window_size
        ub += 1
    end
    while ub > 1 && λs[ub] > λ₀ + window_size
        ub -= 1
    end
    lb, ub
end

"""
    sigma_line(wl)

The cross-section (divided by gf) at wavelength `wl` in Ångstroms of a transition for which the product of the
degeneracy and oscillator strength is `10^log_gf`.
"""
function sigma_line(λ) where F <: Real
    #work in cgs
    e  = electron_charge_cgs
    mₑ = electron_mass_cgs
    c  = c_cgs

    (π*e^2/mₑ/c) * (λ^2/c)
end

"""
    line_profile(λ₀, invΔλ_D, Δλ_L, line_amplitude, λ)

A normalized voigt profile centered on λ₀ with doppler width `invΔλ_D` = 1/Δλ_D and lorentz width 
`Δλ_L` evaluated at `λ` (cm).  Note that this returns values in units of cm^-1.
"""
function line_profile(λ₀, invΔλ_D::F, Δλ_L, line_amplitude::F, λ::F) where F <: Real
    _line_profile(λ₀, invΔλ_D, Δλ_L*invΔλ_D/(4π), line_amplitude*invΔλ_D/sqrt(π), λ)
end
function _line_profile(λ₀, invΔλ_D::F, Δλ_L_invΔλ_D_div_4π::F, amplitude_invΔλ_D_div_sqrt_π::F, λ::F
                      ) where F <: Real
    voigt(Δλ_L_invΔλ_D_div_4π, abs(λ-λ₀) * invΔλ_D) * amplitude_invΔλ_D_div_sqrt_π
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

"""
    voigt(α, v)

The [voigt function](https://en.wikipedia.org/wiki/Voigt_profile#Voigt_functions), ``H``.
"""
function voigt(α, v)
    if α <= 0.2 
        if (v >= 5)
            invv2 = (1/v)^2
            (α/sqrt(π) * invv2) * (1 + 1.5invv2 + 3.75*invv2^2)
        else
            H₀, H₁, H₂ = harris_series(v)
            H₀ + H₁*α + H₂*α^2
        end
    else
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
            r2 = (v*v)/(α*α)
            α_invu = 1/sqrt(2) / ((r2 + 1) * α)
            α2_invu2 = α_invu * α_invu
            sqrt(2/π) * α_invu * (1 + (3*r2 - 1 + ((r2-2)*15*r2+2)*α2_invu2)*α2_invu2)
        end
    end
end
