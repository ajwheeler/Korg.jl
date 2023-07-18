using Interpolations: LinearInterpolation, Flat, lbounds, ubounds
using SpecialFunctions: gamma
using HDF5

"""
    line_absorption(linelist, λs, temp, nₑ, n_densities, partition_fns, ξ
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
- `cuttoff_threshold` (optional, default: 1e-3): see `α_cntm`
"""
function line_absorption!(α, linelist, λs, temp, nₑ, n_densities, partition_fns, ξ, 
                          α_cntm; cutoff_threshold=1e-3)

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
                      n_densities[species"H_I"] .* scaled_vdW.(Ref(line.vdW), m, temp))
        end
        # calculate the lorentz broadening parameter in wavelength. Doing this involves an 
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

        @inbounds view(α, :, lb:ub) .+= line_profile.(line.wl, σ, γ, amplitude, 
                                                      view(concatenated_λs, lb:ub)')
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
at wavelengths `λs` (a vector of ranges).

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
- `use_MHD`: whether or not to use the Mihalas-Daeppen-Hummer formalism to adjust the occupation 
   probabilities of each hydrogen orbital for plasma effects.  Default: `true`.
"""
function hydrogen_line_absorption!(αs, wl_ranges, T, nₑ, nH_I, nHe_I, UH_I, ξ, window_size; 
                                   stark_profiles=_hline_stark_profiles, use_MHD=true)
    λs = vcat(collect.(wl_ranges)...)
    νs = c_cgs ./ λs
    dνdλ = c_cgs ./ λs.^2
    Hmass = get_mass(Formula("H"))

    n_max = maximum(Korg._hline_stark_profiles) do line
        line.upper
    end

    # precalculate occupation probabilities
    ws = if use_MHD
        map(1:n_max) do n
            hummer_mihalas_w(T, n, nH_I, nHe_I, nₑ)
        end
    else 
        ones(n_max)
    end

    β = 1/(kboltz_eV * T)

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
        
        # factor of w because the transition can't happen if the upper level doesn't exist
        levels_factor = ws[line.upper] * (exp(-β*Elo) - exp(-β*Eup)) / UH_I
        amplitude = 10.0^line.log_gf * nH_I * sigma_line(λ₀) * levels_factor

        lb, ub = move_bounds(wl_ranges, 0, 0, λ₀, window_size)
        if lb >= ub
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
    # now do the Brackett series
    n = 4
    for m in 5:n_max
        Elo = Korg.RydbergH_eV * (1 - 1/n^2)
        Eup = Korg.RydbergH_eV * (1 - 1/m^2)
        E = Eup - Elo
        λ0 = Korg.hplanck_eV * Korg.c_cgs / E # cm

        levels_factor = ws[m] * (exp(-β*Elo) - exp(-β*Eup)) / UH_I
        gf = 2 * n^2 * hydrogen_oscillator_strength(n, m)
        amplitude = gf * nH_I * sigma_line(λ0) * levels_factor

        lb, ub = move_bounds(wl_ranges, 0, 0, λ0, window_size)
        lb >= ub && continue
        prof = brackett_line_profile(n, m, view(λs, lb:ub), λ0, T, nₑ, ξ)
        #prof ./= sum(prof) * (λs[lb+1] - λs[lb]) #normalize TODO
        @inbounds view(αs, lb:ub) .+= prof .* amplitude
    end
end

"""
    hydrogen_oscillator_strength(n, m)

The oscillator strength of the transition from the nth to the mth energy level of hydrogen. 
Adapted from HLINOP.f by Peterson and Kurucz. Used for Brackett lines only.

Comparison to the values in [Goldwire 1968](doi.org/10.1086/190180) indicates that this is accurate 
to 10^-4 for the Brackett series.
"""
function hydrogen_oscillator_strength(n, m)
    @assert n < m
    GINF = 0.2027/n^0.71
    GCA = 0.124/n
    FKN = 1.9603n
    WTC = 0.45 - 2.4/n^3 * (n-1)
    FK = FKN*(m/((m-n)*(m+n)))^3
    XMN12 = (m-n)^1.2
    WT = (XMN12 - 1)/(XMN12+WTC)
    FK*(1 - WT*GINF - (0.222+GCA/m)*(1-WT))
end

"""
Normalize stark-broadened line profile (specialized to Brackett series).  Translated and heavily
adapted from HLINOP.f by Peterson and Kurucz via Barklem.

Mostly follows Griem 1960, ApJ, 132, 883, and Griem 1967 with corrections to approximate the Vidal, Cooper & Smith 
(1973, ApJS 25, 37) profiles.

1967 doi.org/10.1086/149097
"""
function brackett_line_profile(n,m,λs,λ₀,T,nₑ,ξ)
    νs = c_cgs ./ λs
    ν₀ = c_cgs / λ₀
    
    ne_1_6 = nₑ^(1/6)
    shielding_parameter = ne_1_6*0.08989/sqrt(T) # the shielding parameter. Called PP in Kurucz
    F0 = 1.25e-9 * nₑ^(2/3) # the Holtzmark field
    GCON1 = 0.2+0.09*sqrt(T/10_0004)/(1+nₑ/1.E13)
    GCON2 = 0.2/(1+nₑ/1.E15)

    # Knm constants as defined by Griem (1960, ApJ 132, 883) for the long 
    # range Holtsmark profile (due to ions only). Lyman and Balmer series 
    # are from VCS, higher series from elsewhere.
    Knm = if (m-n <= 3) && (n<=4)
        Knm_table = [0.0001716 0.0090190 0.1001000 0.5820000
                     0.0005235 0.0177200 0.1710000 0.8660000
                     0.0008912 0.0250700 0.2230000 1.0200000]
        Knm_table[m-n, n]
    else
        # Greim 1960 equation 33 (a=n, b=m, Z=1)
        # 1 / (1 + 0.13(m-n)) is probably a Kurucz addition.  In HLINOP, comment speculates that 
        # it was added to better match the tables.
        5.5e-5 * n^4 * m^4 /(m^2 - n^2) / (1+0.13/(m - n))
    end

    Y1WHT = if (m - n <= 2) && (n <= 2)
        Y1WTM = [1.e18 1e17
                 1.e16 1e14]
        Y1WTM[m-n, n]
    elseif m-n <= 3
       1e14
    else
       1e13
    end
    WTY1 = 1/(1+nₑ/Y1WHT)
    Y1B = 2/(1+0.012/T*sqrt(nₑ/T))
    C1CON = Knm / λ₀ * (m^2 - n^2)^2 / (n^2 * m^2) * 1e-8 # convert factors of cm to Å
    Y1NUM = if m == 2
           550
        elseif (m==3) 
           380
        else
           320
        end
    Y1SCAL = Y1NUM * ((T/10_000)^0.3 / ne_1_6) * WTY1 + Y1B * (1-WTY1)
    C1 = F0*78940/T * C1CON * Y1SCAL

    C2 = F0^2 / (5.96E-23*nₑ) * (Knm/λ₀)^2 * 1e-16 # convert factors of cm to Å

    # Griem 1960 eqn 23.  This is the argument of the Holtsmark profile.
    # β = a / Kₙₘ 
    βs = @. abs(λs-λ₀)/F0/Knm * 1e8 # convert factor of cm to Å

    # y1 is the velocity where the minimum impact parameter and the Lewis cutoff are equal. 
    #   - second order perturbation theory breaks down at the minimum impact parameter
    #   - the impact approximation breaks down at the Lewis cutoff
    # y2 is the where the Lewis cutoff is equal to the Debye length
    # Greim 1967 EQs 6 and 7
    y1 = @. C1*βs
    y2 = @. C2*βs^2

    # TODO: how do these relate to Eqns 6 and 7 in Greim 1967?
    #my_y1 = @. n^2 * hplanck_eV^2 * abs(νs - ν₀) / (2 * kboltz_eV * T)

    G1 = 6.77*sqrt(C1)
    # called F in Kurucz
    impact_electron_contribution = map(y1, y2, βs) do y1, y2, β
        # half-width of the electron impact profile
        # called GAM in Kurucz.  See the  equation between Griem 1967 eqns 13a and 13b.
        hw = if (y2 <= 1e-4) && (y1 <= 1e-5)
            G1*max(0, 0.2114 + log(sqrt(C2)/C1)) * (1-GCON1-GCON2)
        else
            GAM = (G1*(0.5*exp(-min(80, y1))+exponential_integral_1(y1)-0.5*exponential_integral_1(y2))*
                           (1-GCON1/(1+(90*y1)^3)-GCON2/(1+2000*y1)))
            GAM <= 1e-20 ? 0.0 : GAM
        end
        if hw > 0
            # Kurucz source was
            #          F = GAM/PI/(GAM*GAM+BETA*BETA)
            # is this a typo?
            hw/π/(hw*hw+β*β)
        else
            0.0
        end
    end
    
    # called PRQS in Kurucz
    quasistatic_ion_contribution = holtsmark_profile.(βs,shielding_parameter)
    # According to Barklem, a fit eqn 8 of (2nd term) of Griem (1967, ApJ 147, 1092).
    quasistatic_electron_fraction = @. (0.9*y1)^2
    # called FNS in Kurucz
    relative_quasistatic_electron_contribution = (@. (quasistatic_electron_fraction+0.03*sqrt(y1)) / (quasistatic_electron_fraction+1))

    # sqrt(λ/λ₀) corrects the long range part to Δν^(5/2)
    # asymptote, (see Stehle and Hutcheon 1999, A&AS 140, 93).
    dβ_dν = @. c_cgs / (νs^2 * Knm * F0)
    # called STARK1 in Kurucz
    profile = (@. (quasistatic_ion_contribution * (1+relative_quasistatic_electron_contribution) 
                 + impact_electron_contribution)*dβ_dν * sqrt(λs/λ₀))

    # The red wing is multiplied by the Boltzmann factor to roughly account
    # for quantum effects (Stehle 1994, A&AS 104, 509 eqn 7). Assume 
    # absorption case. If emission do for Δν > 0.
    @assert issorted(νs, rev=true)
    i = findlast(νs .< ν₀)
    if !isnothing(i)
        @. profile[begin:i] *= exp((hplanck_cgs*(νs[begin:i] - ν₀))/kboltz_cgs/T)
    end

    # convert from dν to dλ and from cm^-1 to Å^-1
    # TODO where does the factor of 1/2 come from?
    @. 1e8 * c_cgs / λs^2 * profile * 0.5
end

"""
    exponential_integral_1(x)

Compute the first exponential integral, E1(x).  This is a rough approximation lifted from Kurucz's
VCSE1F.
"""
function exponential_integral_1(x)
    if x < 0
        0.0
    elseif x <= 0.01 
        -log(x)-0.577215+x
    elseif x <= 1.0
        -log(x)-0.57721566+x*(0.99999193+x*(-0.24991055+ x*(0.05519968+x*(-0.00976004+ x*0.00107857))))
    elseif x <= 30. 
        (x*(x+2.334733)+0.25062)/(x*(x+3.330657)+ 1.681534)/x*exp(-x)
    else
        0.0
    end
end

const _holtsmark_PROB7 = [ 
     0.005  0.128  0.260  0.389  0.504 
     0.004  0.109  0.220  0.318  0.389
    -0.007  0.079  0.162  0.222  0.244
    -0.018  0.041  0.089  0.106  0.080
    -0.026 -0.003  0.003 -0.023 -0.086
    -0.025 -0.048 -0.087 -0.148 -0.234
    -0.008 -0.085 -0.165 -0.251 -0.343
     0.018 -0.111 -0.223 -0.321 -0.407
     0.032 -0.130 -0.255 -0.354 -0.431
     0.014 -0.148 -0.269 -0.359 -0.427
    -0.005 -0.140 -0.243 -0.323 -0.386
     0.005 -0.095 -0.178 -0.248 -0.307
    -0.002 -0.068 -0.129 -0.187 -0.241
    -0.007 -0.049 -0.094 -0.139 -0.186
    -0.010 -0.036 -0.067 -0.103 -0.143]
const _holtsmark_C7 = [511.318,  1.532,  4.044, 19.266, 41.812]
const _holtsmark_D7 = [-6.070, -4.528, -8.759,-14.984,-23.956]
const _holtsmark_PP = [0.,.2,.4,.6,.8]
const _holtsmark_β_knots = [1.,1.259,1.585,1.995,2.512,3.162,3.981,
                5.012,6.310,7.943,10.,12.59,15.85,19.95,25.12]
"""
    holtsmark_profile(β, P)    

Calculates the Holtsmark profile for broadening of hydrogen lines by quasistatic charged particles.
Adapted from SOFBET in HLINOP by Peterson and Kurucz. Draws heavilly from Griem 1960, ApJ, 132, 883.
"""
function holtsmark_profile(β,P)

    if β > 500 # Very large B
        return (1.5/sqrt(β) + 27/β^2)/β^2
    end

    # Determine relevant Debye range
    # fortran translation note: INT floors floats, not rounds
    IM = min(Int(floor((5*P)+1)),4)
    IP = IM+1
    WTPP = 5*(P-_holtsmark_PP[IM])
    WTPM = 1-WTPP

    if β <= 25.12
        # Indicies into β_boundaries which bound the value of β
        JP = max(2, findfirst(β .<= _holtsmark_β_knots))
        JM = JP - 1
        #println("Js: ", JP, " ", JM)

        #this is linear interpolation into PROB7 wrt β_knots
        WTBP = (β-_holtsmark_β_knots[JM])/(_holtsmark_β_knots[JP]-_holtsmark_β_knots[JM])
        WTBM = 1 - WTBP
        CBP = _holtsmark_PROB7[JP, IP]*WTPP + _holtsmark_PROB7[JP, IM]*WTPM
        CBM = _holtsmark_PROB7[JM, IP]*WTPP + _holtsmark_PROB7[JM, IM]*WTPM
        CORR = 1 + CBP*WTBP + CBM*WTBM

        # Get approximate profile for the inner part
        WT = max(min(0.5*(10-β), 1), 0)

        PR1 = if β <= 10
            PR1 = 8 / (83+(2+0.95*β^2)*β)
        else
            0.0
        end
        PR2 = if β >= 8
            PR2 = (1.5/sqrt(β)+27/β^2)/β^2
        else 
            0.0
        end

        (PR1*WT + PR2*(1 - WT))*CORR
    else
        # Asymptotic part for medium B's (25.12 < B < 500)
        CC=_holtsmark_C7[IP]*WTPP + _holtsmark_C7[IM]*WTPM
        DD=_holtsmark_D7[IP]*WTPP + _holtsmark_D7[IM]*WTPM
        CORR = 1 + DD/(CC + β*sqrt(β))
        (1.5/sqrt(β) + 27/β^2)/β^2*CORR
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
