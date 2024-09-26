using SpecialFunctions: gamma
using ProgressMeter: @showprogress

"""
    line_absorption!(α, linelist, λs, temp, nₑ, n_densities, partition_fns, ξ
                   ; α_cntm=nothing, cutoff_threshold=1e-3, window_size=20.0*1e-8)

Calculate the opacity coefficient, α, in units of cm^-1 from all lines in `linelist`, at wavelengths
`λs` [cm^-1].

other arguments:

  - `temp` the temerature in K (as a vector, for multiple layers, if you like)
  - `n_densities`, a Dict mapping species to absolute number density in cm^-3 (as a vector, if temp is
    a vector).
  - `partition_fns`, a Dict containing the partition function of each species
  - `ξ` is the microturbulent velocity in cm/s (n.b. NOT km/s)
  - `α_cntm` is as a callable returning the continuum opacity as a function of wavelength. The window
    within which a line is calculated will extend to the wavelength at which the Lorentz wings or
    Doppler core of the line are at `cutoff_threshold * α_cntm[line.wl]`, whichever is greater.

# Keyword Arguments

  - `cuttoff_threshold` (default: 3e-4): see `α_cntm`
  - `verbose` (default: false): if true, show a progress bar.
"""
function line_absorption!(α, linelist, λs, temps, nₑ, n_densities, partition_fns, ξ,
                          α_cntm; cutoff_threshold=3e-4, verbose=false)
    if length(linelist) == 0
        return zeros(length(λs))
    end

    # if λs is an vector of ranges, we need this concatenated version for easy indexing
    # the vector of ranges is used for fast index calculations (the move_bounds function).
    concatenated_λs = vcat(λs...)

    #lb and ub are the indices to the upper and lower wavelengths in the "window", i.e. the shortest
    #and longest wavelengths which feel the effect of each line 
    lb = 1
    ub = 1
    β = @. 1 / (kboltz_eV * temps)

    # precompute number density / partition function for each species in the linelist
    n_div_U = map(unique([l.species for l in linelist])) do spec
        spec => @. (n_densities[spec] / partition_fns[spec](log(temps)))
    end |> Dict
    if species"H I" in keys(n_div_U)
        @error "Atomic hydrogen should not be in the linelist. Korg has built-in hydrogen lines."
    end

    # preallocate some arrays for the core loop. 
    # Each element of the arrays corresponds to an atmospheric layer, same at the "temps" array and 
    # the values in "number_densities"
    Γ = Vector{eltype(α)}(undef, size(temps))
    γ = Vector{eltype(α)}(undef, size(temps))
    σ = Vector{eltype(α)}(undef, size(temps))
    amplitude = Vector{eltype(α)}(undef, size(temps))
    levels_factor = Vector{eltype(α)}(undef, size(temps))
    ρ_crit = Vector{eltype(α)}(undef, size(temps))
    inverse_densities = Vector{eltype(α)}(undef, size(temps))
    @showprogress desc="calculating line opacities" enabled=verbose for line in linelist
        m = get_mass(line.species)

        # doppler-broadening width, σ (NOT √[2]σ)
        σ .= doppler_width.(line.wl, temps, m, ξ)

        # sum up the damping parameters.  These are FWHM (γ is usually the Lorentz HWHM) values in 
        # angular, not cyclical frequency (ω, not ν).
        Γ .= line.gamma_rad
        if !ismolecule(line.species)
            @. Γ += nₑ * scaled_stark.(line.gamma_stark, temps)
            Γ .+= n_densities[species"H_I"] .* scaled_vdW.(Ref(line.vdW), m, temps)
        end
        # calculate the lorentz broadening parameter in wavelength. Doing this involves an 
        # implicit aproximation that λ(ν) is linear over the line window.
        # the factor of λ²/c is |dλ/dν|, the factor of 1/2π is for angular vs cyclical freqency,
        # and the last factor of 1/2 is for FWHM vs HWHM
        @. γ = Γ * line.wl^2 / (c_cgs * 4π)

        E_upper = line.E_lower + c_cgs * hplanck_eV / line.wl
        @. levels_factor = exp(-β * line.E_lower) - exp(-β * E_upper)

        #total wl-integrated absorption coefficient
        @. amplitude = 10.0^line.log_gf * sigma_line(line.wl) * levels_factor *
                       n_div_U[line.species]

        ρ_crit .= (line.wl .|> α_cntm) .* cutoff_threshold ./ amplitude
        inverse_densities .= inverse_gaussian_density.(ρ_crit, σ)
        doppler_line_window = maximum(inverse_densities)
        inverse_densities .= inverse_lorentz_density.(ρ_crit, γ)
        lorentz_line_window = maximum(inverse_densities)
        window_size = sqrt(lorentz_line_window^2 + doppler_line_window^2)
        # at present, this line is allocating. Would be good to fix that.
        lb, ub = move_bounds(λs, lb, ub, line.wl, window_size)
        if lb > ub
            continue
        end

        view(α, :, lb:ub) .+= line_profile.(line.wl, σ, γ, amplitude, view(concatenated_λs, lb:ub)')
    end
end

"""
    inverse_gaussian_density(ρ, σ)

Calculate the inverse of a (0-centered) Gaussian PDF with standard deviation `σ`, i.e. the value of
`x` for which `ρ = exp(-0.5 x^2/σ^2}) / √[2π]`, which is given by `σ √[-2 log (√[2π]σρ)]`.  Returns
0 when ρ is larger than any value taken on by the PDF.

See also: [`inverse_lorentz_density`](@ref).
"""
inverse_gaussian_density(ρ, σ) =
    if ρ > 1 / (sqrt(2π) * σ)
        0.0
    else
        σ * sqrt(-2log(sqrt(2π) * σ * ρ))
    end

"""
    inverse_lorentz_density(ρ, γ)

Calculate the inverse of a (0-centered) Lorentz PDF with width `γ`, i.e. the value of `x` for which
`ρ = 1 / (π γ (1 + x^2/γ^2))`, which is given by `√[γ/(πρ) - γ^2]`. Returns 0 when ρ is larger than
any value taken on by the PDF.

See also: [`inverse_gaussian_density`](@ref).
"""
inverse_lorentz_density(ρ, γ) =
    if ρ > 1 / (π * γ)
        0.0
    else
        sqrt(γ / (π * ρ) - γ^2)
    end

"""
    exponential_integral_1(x)

Compute the first exponential integral, E1(x).  This is a rough approximation lifted from Kurucz's
VCSE1F. Used in `brackett_line_profile`.
"""
function exponential_integral_1(x)
    if x < 0
        0.0
    elseif x <= 0.01
        -log(x) - 0.577215 + x
    elseif x <= 1.0
        -log(x) - 0.57721566 +
        x * (0.99999193 + x * (-0.24991055 + x * (0.05519968 + x * (-0.00976004 + x * 0.00107857))))
    elseif x <= 30.0
        (x * (x + 2.334733) + 0.25062) / (x * (x + 3.330657) + 1.681534) / x * exp(-x)
    else
        0.0
    end
end

"""
    doppler_width(λ₀ T, m, ξ)

The standard deviation of of the doppler-broadening profile.  In standard spectroscopy texts, the
Doppler width often refers to σ√2, but this is σ
"""
doppler_width(λ₀, T, m, ξ) = λ₀ * sqrt(kboltz_cgs * T / m + (ξ^2) / 2) / c_cgs

"""
the stark broadening gamma scaled acording to its temperature dependence
"""
scaled_stark(γstark, T; T₀=10_000) = γstark * (T / T₀)^(1 / 6)

"""
    scaled_vdW(vdW, m, T)

The vdW broadening gamma scaled acording to its temperature dependence, using either simple scaling
or ABO. See Anstee & O'Mara (1995) or https://www.astro.uu.se/~barklem/howto.html for the definition
of the ABO γ.

`vdW` should be either `γ_vdW` evaluated at 10,000 K, or tuple containing the ABO params `(σ, α)`.
The species mass, `m`, is ignored in the former case.
"""
scaled_vdW(vdW::Real, m, T, T₀=10_000) = vdW * (T / T₀)^0.3
function scaled_vdW(vdW::Tuple{F,F}, m, T) where F<:Real
    v₀ = 1e6 #σ is given at 10_000 m/s = 10^6 cm/s
    σ = vdW[1]
    α = vdW[2]

    invμ = 1 / (1.008 * amu_cgs) + 1 / m #inverse reduced mass
    vbar = sqrt(8 * kboltz_cgs * T / π * invμ) #relative velocity
    #n.b. "gamma" is the gamma function, not a broadening parameter
    2 * (4 / π)^(α / 2) * gamma((4 - α) / 2) * v₀ * σ * (vbar / v₀)^(1 - α)
end

"""
    sigma_line(wl)

The cross-section (divided by gf) at wavelength `wl` in Ångstroms of a transition for which the
product of the degeneracy and oscillator strength is `10^log_gf`.
"""
function sigma_line(λ::Real)
    #work in cgs
    e = electron_charge_cgs
    mₑ = electron_mass_cgs
    c = c_cgs

    #the factor of |dλ/dν| = λ²/c is because we are working in wavelength rather than frequency
    (π * e^2 / mₑ / c) * (λ^2 / c)
end

"""
    line_profile(λ₀, σ, γ, amplitude, λ)

A voigt profile centered on λ₀ with Doppler width σ (NOT √[2] σ, as the "Doppler width" is often
defined) and Lorentz HWHM γ evaluated at `λ` (cm).  Returns values in units of cm^-1.
"""
function line_profile(λ₀::Real, σ::Real, γ::Real, amplitude::Real, λ::Real)
    inv_σsqrt2 = 1 / (σ * sqrt(2))
    scaling = inv_σsqrt2 / sqrt(π) * amplitude
    voigt_hjerting(γ * inv_σsqrt2, abs(λ - λ₀) * inv_σsqrt2) * scaling
end

@inline function harris_series(v) # assume v < 5
    v2 = v * v
    H₀ = exp(-(v2))
    H₁ = if (v < 1.3)
        -1.12470432 + (-0.15516677 + (3.288675912 + (-2.34357915 + 0.42139162 * v) * v) * v) * v
    elseif v < 2.4
        -4.48480194 + (9.39456063 + (-6.61487486 + (1.98919585 - 0.22041650 * v) * v) * v) * v
    else #v < 5
        ((0.554153432 +
          (0.278711796 + (-0.1883256872 + (0.042991293 - 0.003278278 * v) * v) * v) * v) /
         (v2 - 3 / 2))
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
    v2 = v * v
    if α <= 0.2 && (v >= 5)
        invv2 = (1 / v2)
        (α / sqrt(π) * invv2) * (1 + 1.5invv2 + 3.75 * invv2^2)
    elseif α <= 0.2 #v < 5
        H₀, H₁, H₂ = harris_series(v)
        H₀ + (H₁ + H₂ * α) * α
    elseif (α <= 1.4) && (α + v < 3.2)
        #modified harris series: M_i is H'_i in the source text
        H₀, H₁, H₂ = harris_series(v)
        M₀ = H₀
        M₁ = H₁ + 2 / sqrt(π) * M₀
        M₂ = H₂ - M₀ + 2 / sqrt(π) * M₁
        M₃ = 2 / (3sqrt(π)) * (1 - H₂) - (2 / 3) * v2 * M₁ + (2 / sqrt(π)) * M₂
        M₄ = 2 / 3 * v2 * v2 * M₀ - 2 / (3sqrt(π)) * M₁ + 2 / sqrt(π) * M₃
        ψ = 0.979895023 + (-0.962846325 + (0.532770573 - 0.122727278 * α) * α) * α
        ψ * (M₀ + (M₁ + (M₂ + (M₃ + M₄ * α) * α) * α) * α)
    else #α > 1.4 or (α > 0.2 and α + v > 3.2)
        r2 = (v2) / (α * α)
        α_invu = 1 / sqrt(2) / ((r2 + 1) * α)
        α2_invu2 = α_invu * α_invu
        sqrt(2 / π) * α_invu * (1 + (3 * r2 - 1 + ((r2 - 2) * 15 * r2 + 2) * α2_invu2) * α2_invu2)
    end
end
