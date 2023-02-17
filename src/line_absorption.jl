using Interpolations: LinearInterpolation, Flat, lbounds, ubounds
using SpecialFunctions: gamma
using HDF5

"""
    line_absorption(linelist, λs, temp, nₑ, n_densities, partition_fns, ξ
                   ; α_cntm=nothing, cutoff_threshold=1e-3, window_size=20.0*1e-8)

Calculate the opacity coefficient, α, in units of cm^-1 from all lines in `linelist`, at wavelengths
`λs` [cm^-1]. 

other arguments:
- `temp` the temerature in K (at multiply layers, if you like)
- `n_densities`, a Dict mapping species to absolute number density [cm^-3] (as a vector, if temp is
   a vector).
- `partition_fns`, a Dict containing the partition function of each species
- `ξ` is the microturbulent velocity in cm/s (n.b. NOT km/s)
- `α_cntm` is as a callable returning the continuum opacity as a function of wavelength. The window 
   within which a line is calculated will extend to the wavelength at which the Lorentz wings or 
   Doppler core of the line are at `cutoff_threshold * α_cntm[line.wl]`, whichever is greater.  
- `cuttoff_threshold` (optional, default: 1e-3): see `α_cntm`
"""
function line_absorption!(α, linelist, λs, temp, nₑ, n_densities, partition_fns, ξ, 
                          α_cntm; cutoff_threshold=1e-3)
    if length(linelist) == 0
        return zeros(length(λs))
    end

    #lb and ub are the indices to the upper and lower wavelengths in the "window", i.e. the shortest
    #and longest wavelengths which feel the effect of each line 
    lb = 1
    ub = 1
    β = @. 1/(kboltz_eV * temp)

    # precompute number density / partition function for each species in the linelist
    n_div_Z = map(collect(Set([l.species for l in linelist]))) do spec
        spec => @. (n_densities[spec] / partition_fns[spec](log(temp)))
    end |> Dict

    for line in linelist
        m = get_mass(line.species)
        
        # doppler-broadening width, σ (NOT √[2]σ)
        σ = doppler_width.(line.wl, temp, m, ξ)

        # sum up the damping parameters.  These are FWHM (γ is usually the Lorentz HWHM) values in 
        # angular, not cyclical frequency (ω, not ν).
        Γ = line.gamma_rad 
        if !ismolecule(line.species) 
            Γ = Γ .+ (nₑ .* scaled_stark.(line.gamma_stark, temp) +
                      n_densities[species"H_I"] .* [scaled_vdW(line.vdW, m, T) for T in temp])
        end
        # calculate the lorentz broadening parameter in in wavelength. Doing this involves an 
        # implicit aproximation that λ(ν) is linear over the line window.
        # the factor of λ²/c is |dλ/dν|, the factor of 1/2π is for angular vs cyclical freqency,
        # and the last factor of 1/2 is for FWHM vs HWHM
        γ = @. Γ * line.wl^2 / (c_cgs * 4π)

        E_upper = line.E_lower + c_cgs * hplanck_eV / line.wl 
        levels_factor = (@. (exp(-β*line.E_lower) - exp(-β*E_upper)))

        #total wl-integrated absorption coefficient
        amplitude = @. 10.0^line.log_gf*sigma_line(line.wl)*levels_factor*n_div_Z[line.species]

        ρ_crit = [cntm(line.wl) * cutoff_threshold for cntm in α_cntm] ./ amplitude
        doppler_line_window = maximum(inverse_gaussian_density.(ρ_crit, σ))
        lorentz_line_window = maximum(inverse_lorentz_density.(ρ_crit, γ))
        window_size = sqrt(lorentz_line_window^2 + doppler_line_window^2)
        lb, ub = move_bounds(λs, lb, ub, line.wl, window_size)
        if lb > ub
            continue
        end

        @inbounds view(α, :, lb:ub) .+= line_profile.(line.wl, σ, γ, amplitude, view(λs, lb:ub)')
    end
end

"""
    inverse_gaussian_density(ρ, σ)

Calculate the inverse of a (0-centered) Gaussian PDF with standard deviation `σ`, i.e. the value of 
`x` for which `ρ = exp(-0.5 x^2/σ^2}) / √[2π]`, which is given by `σ √[-2 log (√[2π]σρ)]`.  Returns 
0 when ρ is larger than any value taken on by the PDF.

See also: [`inverse_lorentz_density`](@ref).
"""
function inverse_gaussian_density(ρ, σ)
    if ρ > 1/(sqrt(2π) * σ)
        0.0
    else
        σ * sqrt(-2log(sqrt(2π) * σ * ρ))
    end
end

"""
    inverse_lorentz_density(ρ, γ)

Calculate the inverse of a (0-centered) Lorentz PDF with width `γ`, i.e. the value of `x` for which 
`ρ = 1 / (π γ (1 + x^2/γ^2))`, which is given by `√[γ/(πρ) - γ^2]`. Returns 0 when ρ is larger than 
any value taken on by the PDF.

See also: [`inverse_gaussian_density`](@ref).
"""
function inverse_lorentz_density(ρ, γ)
    if ρ > 1/(π*γ)
        0.0
    else
        sqrt(γ/(π * ρ) - γ^2)
    end
end

#load Stark broadening profiles from disk
function _load_stark_profiles(fname)
    h5open(fname, "r") do fid
        map(keys(fid)) do transition
            temps = read(fid[transition], "temps")
            nes = read(fid[transition], "electron_number_densities")
            delta_nu_over_F0 = read(fid[transition], "delta_nu_over_F0")
            P = read(fid[transition], "profile")

            logP = log.(P)
            # clipping these to finite numbers avoid nans when interpolating
            # -700 is slightly larger than log(-floatmax())
            logP[logP .== -Inf] .= -700 
            
            (temps=temps, 
             electron_number_densities=nes,
             #flat BCs to allow wls very close to the line center
             profile=LinearInterpolation((temps, nes, [-floatmax() ; log.(delta_nu_over_F0[2:end])]), 
                                         logP; extrapolation_bc=Flat()),
             lower = HDF5.read_attribute(fid[transition], "lower"),
             upper = HDF5.read_attribute(fid[transition], "upper"),
             Kalpha = HDF5.read_attribute(fid[transition], "Kalpha"),
             log_gf = HDF5.read_attribute(fid[transition], "log_gf"),
             λ0 = LinearInterpolation((temps, nes), read(fid[transition], "lambda0") * 1e-8) 
            )
        end
    end
end
const _hline_stark_profiles = _load_stark_profiles(joinpath(_data_dir, 
                                                      "Stehle-Hutchson-hydrogen-profiles.h5"))

#used in hydrogen_line_absorption
_zero2epsilon(x) = x + (x == 0) * floatmin()

"""
    hydrogen_line_absorption!(αs, λs, T, nₑ, nH_I, UH_I, ξ, window_size; kwargs...)

Calculate contribution to the the absorption coefficient, αs, from hydrogen lines in units of cm^-1,
at wavelengths `λs`.  TODO MHD

Uses profiles from [Stehlé & Hutcheon (1999)](https://ui.adsabs.harvard.edu/abs/1999A%26AS..140...93S/abstract),
which include Stark and Doppler broadening.  
For Halpha, Hbeta, and Hgamma, the p-d approximated profiles from 
[Barklem, Piskunovet, and O'Mara 2000](https://ui.adsabs.harvard.edu/abs/2000A%26A...363.1091B/abstract)
are added to the absortion coefficient.  This "convolution by summation" is inexact, but true 
convolution is expensive.

Arguments:
- `T`: temperature [K]
- `nₑ`: electron number density [cm^-3]
- `nH_I`: neutral hydrogen number density [cm^-3]
- `UH_I`: the value of the neutral hydrogen partition function
- `ξ`: microturbulent velocity [cm/s]. This is only applied to Hα-Hγ.  Other hydrogen lines profiles 
   are dominated by stark broadening, and the stark broadened profiles are pre-convolved with a 
   doppler profile.
- `window_size`: the max distance from each line center [cm] at which to calculate the stark
   and self broadening profiles
   absorption for Hα-Hγ (those dominated by self-broadening).

Keyword arguments:
- `stark_profiles` (default: `Korg._hline_stark_profiles`): tables from which to interpolate Stark 
   profiles
- `use_MHD` TODO
"""
function hydrogen_line_absorption!(αs, λs, T, nₑ, nH_I, nHe_I, UH_I, ξ, window_size; 
                                   stark_profiles=_hline_stark_profiles, use_MHD=true)
    νs = c_cgs ./ λs
    dνdλ = c_cgs ./ λs.^2
    Hmass = get_mass(Formula("H"))

    n_max = maximum(Korg._hline_stark_profiles) do line
        line.upper
    end

    # precalculate occupation probabilities
    ws = map(1:n_max) do n
        E = RydbergH_eV * (1 - 1/n^2)
        hummer_mihalas_w(T, n, nH_I, nHe_I, nₑ)
    end

    #This is the Holtzmark field, by which the frequency-unit-detunings are divided for the 
    #interpolated stark profiles
    F0 = 1.25e-9 * nₑ^(2/3)
    for line in stark_profiles
        if !all(lbounds(line.λ0.itp)[1:2] .< (T, nₑ) .< ubounds(line.λ0.itp)[1:2])
            continue #transitions to high levels are omitted for high nₑ and T
        end
        λ₀ = line.λ0(T, nₑ)

        Elo = RydbergH_eV * (1 - 1/line.lower^2)
        Eup = RydbergH_eV * (1 - 1/line.upper^2)
        β = 1/(kboltz_eV * T)
        
        levels_factor = if use_MHD
            # the transition can't happen if the upper level doesn't exist
            ws[line.upper] * (exp(-β*Elo) - exp(-β*Eup)) / UH_I
        else
            (exp(-β*Elo) - exp(-β*Eup)) / UH_I
        end
        amplitude = 10.0^line.log_gf * nH_I * sigma_line(λ₀) * levels_factor

        lb, ub = move_bounds(λs, 0, 0, λ₀, window_size)
        if lb == ub
            continue
        end
         
        # if it's Halpha, Hbeta, or Hgamma, add the resonant broadening to the absorption vector
        # use the Barklem+ 2000 p-d approximation
        if line.lower == 2 && line.upper in [3, 4, 5]
            #ABO params and line center
            λ₀, σABO, αABO = if line.upper == 3
                6.56460998e-5, 1180.0, 0.677
            elseif line.upper == 4
                4.8626810200000004e-5, 2320.0, 0.455
            elseif line.upper == 5
                4.34168232e-5, 4208.0, 0.380
            end

            Γ = scaled_vdW((σABO*bohr_radius_cgs^2, αABO), Hmass, T) * nH_I
            # convert to HWHM wavelength units. (see comment in line_absorption! for explanation)
            γ = Γ * λ₀^2 / (c_cgs * 4π) 

            σ = doppler_width(λ₀, T, Hmass, ξ)

            @inbounds view(αs,lb:ub) .+= line_profile.(λ₀, σ, γ, amplitude, view(λs, lb:ub))
        end

        # Stehle+ 1999 Stark-broadened profiles
        ν₀ = c_cgs / (λ₀)
        scaled_Δν = _zero2epsilon.(abs.(view(νs,lb:ub) .- ν₀) ./ F0)
        dIdν = exp.(line.profile.(T, nₑ, log.(scaled_Δν)))
        @inbounds view(αs, lb:ub) .+= dIdν .* view(dνdλ, lb:ub) .* amplitude
    end
end

"""
    doppler_width(λ₀ T, m, ξ)

The standard deviation of of the doppler-broadening profile.  In standard spectroscopy texts, the 
Doppler width often refers to σ√2, but this is σ
"""
doppler_width(λ₀, T, m, ξ) = λ₀ * sqrt(kboltz_cgs*T / m + (ξ^2)/2) / c_cgs

"the stark broadening gamma scaled acording to its temperature dependence"
scaled_stark(γstark, T; T₀=10_000) = γstark * (T/T₀)^(1/6)
             
"""
    scaled_vdW(vdW, m, T)

The vdW broadening gamma scaled acording to its temperature dependence, using either simple scaling 
or ABO. See Anstee & O'Mara (1995) or https://www.astro.uu.se/~barklem/howto.html for the definition
of the ABO γ. 

`vdW` should be either `γ_vdW` evaluated at 10,000 K, or tuple containing the ABO params `(σ, α)`. 
The species mass, `m`, is ignored in the former case.
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

"""
    sigma_line(wl)

The cross-section (divided by gf) at wavelength `wl` in Ångstroms of a transition for which the 
product of the degeneracy and oscillator strength is `10^log_gf`.
"""
function sigma_line(λ::Real)
    #work in cgs
    e  = electron_charge_cgs
    mₑ = electron_mass_cgs
    c  = c_cgs

    #the factor of |dλ/dν| = λ²/c is because we are working in wavelength rather than frequency
    (π*e^2/mₑ/c) * (λ^2/c)
end

"""
    line_profile(λ₀, σ, γ, amplitude, λ)

A voigt profile centered on λ₀ with Doppler width σ (NOT √[2] σ, as the "Doppler width" is often 
defined) and Lorentz HWHM γ evaluated at `λ` (cm).  Returns values in units of cm^-1.
"""
function line_profile(λ₀::Real, σ::Real, γ::Real, amplitude::Real, λ::Real)
    inv_σsqrt2 = 1/(σ*sqrt(2))
    scaling = inv_σsqrt2 / sqrt(π) * amplitude
    voigt_hjerting(γ*inv_σsqrt2, abs(λ-λ₀)*inv_σsqrt2) * scaling
end

@inline function harris_series(v) # assume v < 5
    v2 = v*v
    H₀ = exp(-(v2))
    H₁ = if (v < 1.3)
        -1.12470432 + (-0.15516677 + (3.288675912 +(-2.34357915 + 0.42139162*v)*v)*v)*v
    elseif v < 2.4
        -4.48480194 + (9.39456063 + (-6.61487486 + (1.98919585  - 0.22041650*v)*v)*v)*v
    else #v < 5
        ((0.554153432 + (0.278711796 + (-0.1883256872 + (0.042991293 - 0.003278278*v)*v)*v)*v) / 
         (v2 - 3/2))
    end
    H₂ = (1 - 2v2) * H₀
    H₀, H₁, H₂
end

"""
    voigt_hjerting(α, v)

The [Hjerting function](https://en.wikipedia.org/wiki/Voigt_profile#Voigt_functions), ``H``, 
somtimes called the Voigt-Hjerting function. ``H`` is defined as
`H(α, v) = ∫^∞_∞ exp(-y^2) / ((u-y)^2 + α^2) dy`
(see e.g. the unnumbered equation after Gray equation 11.47).  It is equal to the ratio of the 
absorption coefficient to the value of the absorption coefficient obtained at the line center with 
only Doppler broadening.

If `x = λ-λ₀`, `Δλ_D = σ√2` is the Doppler width, and `Δλ_L = 4πγ` is the Lorentz width,
```
voigt(x|Δλ_D, Δλ_L) = H(Δλ_L/(4πΔλ_D), x/Δλ_D) / (Δλ_D√π)
                    = H(γ/(σ√2), x/(σ√2)) / (σ√(2π))
```
Approximation from [Hunger 1965](https://ui.adsabs.harvard.edu/abs/1956ZA.....39...36H/abstract).
"""
function voigt_hjerting(α, v)
    v2 = v*v
    if α <= 0.2 && (v >= 5)
        invv2 = (1/v2)
        (α/sqrt(π) * invv2) * (1 + 1.5invv2 + 3.75*invv2^2)
    elseif α <= 0.2 #v < 5
        H₀, H₁, H₂ = harris_series(v)
        H₀ + (H₁ + H₂*α)*α
    elseif (α <= 1.4) && (α + v < 3.2)
        #modified harris series: M_i is H'_i in the source text
        H₀, H₁, H₂ = harris_series(v)
        M₀ = H₀
        M₁ = H₁ + 2/sqrt(π) * M₀
        M₂ = H₂ - M₀ + 2/sqrt(π) * M₁
        M₃ = 2/(3sqrt(π))*(1-H₂) - (2/3)*v2*M₁ + (2/sqrt(π))*M₂
        M₄ = 2/3 * v2*v2 * M₀ - 2/(3sqrt(π)) * M₁ + 2/sqrt(π) * M₃
        ψ =  0.979895023 + (-0.962846325 + (0.532770573 - 0.122727278*α)*α)*α
        ψ * (M₀ + (M₁ + (M₂ + (M₃ + M₄*α)*α)*α)*α)
    else #α > 1.4 or (α > 0.2 and α + v > 3.2)
        r2 = (v2)/(α*α)
        α_invu = 1/sqrt(2) / ((r2 + 1) * α)
        α2_invu2 = α_invu * α_invu
        sqrt(2/π) * α_invu * (1 + (3*r2 - 1 + ((r2-2)*15*r2+2)*α2_invu2)*α2_invu2)
    end
end
