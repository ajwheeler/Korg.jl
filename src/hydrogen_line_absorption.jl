using HDF5
using SpecialFunctions: gamma
# LinearInterpolation is now deprecated in favor of linear_interpolation, but we'll keep using the 
# old version for a while for compatibility
using Interpolations: LinearInterpolation, lbounds, ubounds, Flat

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
            
            profile = LinearInterpolation((temps, nes, [-floatmax() ; log.(delta_nu_over_F0[2:end])]), 
                                          logP;
                                          extrapolation_bc=Flat())

            λ0 = LinearInterpolation((temps, nes), read(fid[transition], "lambda0") * 1e-8) 
             
            (
             temps=temps, 
             electron_number_densities=nes,
             profile=profile,
             #flat BCs to allow wls very close to the line center
             lower = HDF5.read_attribute(fid[transition], "lower"),
             upper = HDF5.read_attribute(fid[transition], "upper"),
             Kalpha = HDF5.read_attribute(fid[transition], "Kalpha"),
             log_gf = HDF5.read_attribute(fid[transition], "log_gf"),
             λ0 = λ0
            )
        end
    end
end
const _hline_stark_profiles = _load_stark_profiles(
    joinpath(_data_dir, "Stehle-Hutchson-hydrogen-profiles.h5"))

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
- `window_size`: the max distance from each line center [cm] at which to calculate the Stark and 
   self-broadening profiles absorption for Hα-Hγ (those dominated by self-broadening).

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

    n_max = maximum(_hline_stark_profiles) do line
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

    #This is the Holtsmark field, by which the frequency-unit-detunings are divided for the 
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
    E_low = RydbergH_eV * (1 - 1/n^2)
    for m in 5:n_max
        E = RydbergH_eV * (1/n^2 - 1/m^2)
        λ0 = hplanck_eV * c_cgs / E # cm
        levels_factor = ws[m] * exp(-β*E_low) * (1 - exp(-β*E)) / UH_I
        gf = 2 * n^2 * brackett_oscillator_strength(n, m)
        amplitude = gf * nH_I * sigma_line(λ0) * levels_factor

        stark_profile_itp, stark_window = bracket_line_interpolator(m, λ0, T, nₑ, ξ, wl_ranges[1][1], wl_ranges[end][end])
        lb, ub = move_bounds(wl_ranges, 0, 0, λ0, stark_window)

        # renormalize profile?
        @inbounds view(αs, lb:ub) .+= stark_profile_itp.(view(λs, lb:ub)) .* amplitude
    end

end

"""
    brackett_oscillator_strength(n, m)

The oscillator strength of the transition from the 4th to the mth energy level of hydrogen. 
Adapted from HLINOP.f by Peterson and Kurucz. 

Comparison to the values in [Goldwire 1968](https://doi.org/10.1086/190180) indicates that this is 
accurate to 10^-4 for the Brackett series.
"""
function brackett_oscillator_strength(n, m)
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
    bracket_line_interpolator(m, λ₀, T, nₑ, ξ, λmin, λmax; kwargs...)

This routine numerically convolves the two components of the Brackett line stark profile 
(quasistatic/Holtsmark and impact) and the doppler profile, if necessary.  It returns a pair 
containing the interpolator and the distance from the line center at which it is defined.

# Arguments
- `m`: the principle quantum number of the upper level
- `λ₀`: the line center in Å
- `T`: the temperature [K]
- `nₑ`: the electron number density [cm^-3]
- `ξ`: the microturbulence [cm/s]
- `λmin`: the minimum wavelength at which the profile should be computed (used to avoid calculations 
   outside the required region)
- `λmax`: the mxinimum wavelength at which the profile should be computed 

# Keyword Arguments
- `n_wavelength_points` (default=201): the number of wavelengths at which to sample the profiles
  quasistatic profiles to be convolved.
- `window_size` (default=5): the size of the wavelength range over which the profiles should be 
 calculated, in units of the characteristic profile width
"""
function bracket_line_interpolator(m, λ₀, T, nₑ, ξ, λmin=0, λmax=Inf; 
                                   n_wavelength_points=201, window_size=5, 
                                   include_doppler_threshold=0.25)
    n = 4

    # get stark width
    F0 = 1.25e-9 * nₑ^(2/3) # the Holtsmark field
    Knm = greim_1960_Knm(n, m)
    stark_width = 1.6678E-18 * Knm * F0 * c_cgs

    # get doppler width
    σdop = doppler_width(λ₀, T, atomic_masses[1], ξ)

    # set wavelengths for calculations and convolutions
    window = window_size * max(σdop, stark_width)
    λstart = max(λmin, λ₀ - window)
    λend = min(λ₀ + window, λmax)
    if λstart > λmax || λend < λmin || λstart == λend # 3rd case for one-wavelength synthesis
        # if the calculated wavelength window is entirely outside the synthesis range, return a noop
        # interpolator and a null window (for type stability)
        return LinearInterpolation([], []), 0.0
    end
    wls = range(λstart, λend; length=n_wavelength_points)
    start_ind = (n_wavelength_points-1) ÷ 2 # used to get indices corresponding to original wls

    # compute stark profiles
    ϕ_impact, ϕ_quasistatic = brackett_line_stark_profiles(m, wls, λ₀, T, nₑ)

    # possibly include doppler by convolving it with the quasistatic profile
    if σdop / stark_width > include_doppler_threshold
        ϕ_dop = normal_pdf.(wls .- λ₀, σdop)
        ϕ_quasistatic = autodiffable_conv(ϕ_quasistatic, ϕ_dop) * step(wls)
        ϕ_quasistatic = ϕ_quasistatic[start_ind:start_ind+n_wavelength_points-1]
    end

    # convolve impact and quasistatic profiles
    ϕ_conv = autodiffable_conv(ϕ_impact, ϕ_quasistatic) * step(wls)

    itp = LinearInterpolation(wls, ϕ_conv[start_ind:start_ind+n_wavelength_points-1])
    itp, window
end

"""
    brackett_line_stark_profiles(m, λs, λ₀, T, nₑ)

Stark-broadened line profile (specialized to Brackett series).  Translated and heavily
adapted from HLINOP.f by Barklem, who adapted it from Peterson and Kurucz.  Mostly follows 
[Griem 1960](https://doi.org/10.1086/146987), and [Griem 1967](https://doi.org/10.1086/149097).  Ions and distant
electrons have E fields which can be treated quasi-statically, leading to a
[Holtsmark broadening profile](https://doi.org/10.1002/andp.19193630702).

Returns a pair of vectors containing the impact and quasistatic profiles, respectively.

Arguents:
- `m`: the upper level of the transition
- `λs`: the wavelengths at which to calculate the profile [cm]
- `λ₀`: the line center [cm]
- `T`: temperature [K]
- `nₑ`: electron number density [cm^-3]
"""
function brackett_line_stark_profiles(m, λs, λ₀, T, nₑ)
    n = 4 # Brackett lines only
    νs = c_cgs ./ λs
    ν₀ = c_cgs / λ₀
    
    ne_1_6 = nₑ^(1/6)
    F0 = 1.25e-9 * nₑ^(2/3) # the Holtsmark field
    GCON1 = 0.2+0.09*sqrt(T/10_0004)/(1+nₑ/1.E13)
    GCON2 = 0.2/(1+nₑ/1.E15)

    Knm = greim_1960_Knm(n, m)

    Y1WHT = if m-n <= 3 
       1e14
    else
       1e13
    end
    WTY1 = 1/(1+nₑ/Y1WHT)
    Y1B = 2/(1+0.012/T*sqrt(nₑ/T))
    C1CON = Knm / λ₀ * (m^2 - n^2)^2 / (n^2 * m^2) * 1e-8 # convert factors of cm to Å
    Y1NUM = 320 # specialized to m=4
    Y1SCAL = Y1NUM * ((T/10_000)^0.3 / ne_1_6) * WTY1 + Y1B * (1-WTY1)
    C1 = F0*78940/T * C1CON * Y1SCAL

    C2 = F0^2 / (5.96E-23*nₑ) * (Knm/λ₀)^2 * 1e-16 # convert factors of cm to Å

    # Griem 1960 eqn 23.  This is the argument of the Holtsmark profile.
    βs = @. abs(λs-λ₀)/F0/Knm * 1e8 # convert factor of cm to Å, for the calculations below

    # y, introduced in Griem 1967 EQ 5, is the ratio of particle kinetic energy to average kinetic 
    # energy: (1/2 mv^2) / (kT)
    # y1 corresponds to velocity where the minimum impact parameter and the Lewis cutoff are equal. 
    #   - second order perturbation theory breaks down at the minimum impact parameter
    #   - the impact approximation breaks down at the Lewis cutoff
    # y2 is the where the Lewis cutoff is equal to the Debye length
    # Greim 1967 EQs 6 and 7 defines these quantities, but I'm slightly confused about how they
    # are related to these definitions.
    y1 = @. C1*βs
    y2 = @. C2*βs^2

    G1 = 6.77*sqrt(C1)
    # called F in Kurucz
    impact_electron_profile = map(y1, y2, βs) do y1, y2, β
        # width of the electron impact profile
        # called GAM in Kurucz.  See the equation between Griem 1967 eqns 13a and 13b.
        width = if (y2 <= 1e-4) && (y1 <= 1e-5)
            G1*max(0, 0.2114 + log(sqrt(C2)/C1)) * (1-GCON1-GCON2)
        else
            GAM = (G1*(0.5*exp(-min(80, y1))+exponential_integral_1(y1)-0.5*exponential_integral_1(y2))*
                           (1-GCON1/(1+(90*y1)^3)-GCON2/(1+2000*y1)))
            GAM <= 1e-20 ? 0.0 : GAM
        end

        if width > 0
            width / (π*(width^2 + β^2)) # Lorentz density
        else
            zero(typeof(β)) #make it type-stable
        end
    end
    
    shielding_parameter = ne_1_6*0.08989/sqrt(T) # the shielding parameter. Called PP in Kurucz
    quasistatic_ion_contribution = holtsmark_profile.(βs,shielding_parameter) # called PRQS in Kurucz


    # quasistatic_e_contrib is a fit to (sqrt(π) - 2*gamma(3/2, y1))/sqrt(π) (taken from HLINOP/Kurucz), 
    # the second term in eqn 8 of Griem (1967, ApJ 147, 1092). The sum in Greim is an expansion of 
    # the gamma function. 
    ps = (0.9*y1).^2
    quasistatic_e_contrib = @. (ps+0.03*sqrt(y1))/(ps+1.)

    total_quasistatic_profile = @. quasistatic_ion_contribution * (1+quasistatic_e_contrib) 
    # this correction makes the profile not normalized.  It's unclear to me that we should be 
    # scaling the profile and not the number of perturbers used to calculate it.

    dβ_dλ = 1e8 / (Knm * F0)
    # is this the appropriate treatment? document.
    for profile in [impact_electron_profile, total_quasistatic_profile]
        # sqrt(λ/λ₀) corrects the long range part to Δν^(5/2)
        # asymptote, (see Stehle and Hutcheon 1999, A&AS 140, 93).
        @. profile .*= sqrt(λs/λ₀)

        # The red wing is multiplied by the Boltzmann factor to roughly account
        # for quantum effects (Stehle 1994, A&AS 104, 509 eqn 7). Assume 
        # absorption case. If emission do for Δν > 0.
        @assert issorted(νs, rev=true)
        i = findlast(νs .< ν₀)
        if !isnothing(i)
            profile[begin:i] .*= @. exp((hplanck_cgs*(νs[begin:i] - ν₀))/kboltz_cgs/T)
        end

        profile .*= dβ_dλ
    end

    impact_electron_profile, total_quasistatic_profile
end

const _greim_Kmn_table = [
    0.0001716 0.0090190 0.1001000 0.5820000
    0.0005235 0.0177200 0.1710000 0.8660000
    0.0008912 0.0250700 0.2230000 1.0200000
    ]

"""
    greim_1960_Knm(n, m)

Knm constants as defined by [Griem 1960](https://doi.org/10.1086/146987) for the long range Holtsmark 
profile. This function includes only the values for Brackett lines.

``K_{nm} = C_{nm} 2 π c / λ^2`` where ``C_{nm} F = Δω`` and ``F`` is the ion field.  
See Griem 1960 EQs 7 and 12.
This works out to ``K_{nm} = λ/F``.
"""
function greim_1960_Knm(n, m)
    if (m-n <= 3) && (n<=4)
        _greim_Kmn_table[m-n, n]
    else
        # Greim 1960 equation 33 (a=n, b=m, Z=1)
        # 1 / (1 + 0.13(m-n)) is probably a Kurucz addition.  In Barklem's HLINOP, comment 
        # speculates that it was added to better match the tables.
        5.5e-5 * n^4 * m^4 /(m^2 - n^2) / (1+0.13/(m - n))
    end
end

"""
    holtsmark_profile(β, P)    

Calculates the Holtsmark profile for broadening of hydrogen lines by quasistatic charged particles.
Adapted from SOFBET in HLINOP by Peterson and Kurucz. Draws heavily from 
[Griem 1960](https://doi.org/10.1086/146987).
"""
function holtsmark_profile(β,P)

    if β > 500 # Very large β
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

# used in holtsmark_profile
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

