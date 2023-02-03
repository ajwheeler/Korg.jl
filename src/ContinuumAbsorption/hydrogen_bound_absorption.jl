using Interpolations: LinearInterpolation, Line
using HDF5 #TODO remove?

# load hydrogen bf cross sections
_H_cross_sections = let
    _H_level_ns, _H_level_Ls, _H_level_binding_Es, _H_level_cross_sections = 
            h5open(joinpath(_data_dir, "bf_cross-sections", "individual_H_cross-sections.h5")) do f
        sigmas = map(eachcol(read(f["E"])), eachcol(read(f["sigma"]))) do Es, σs
            LinearInterpolation(Es, σs, extrapolation_bc=Line())
        end
        read(f["n"]), read(f["L"]), read(f["E_bind"]), sigmas
    end

    zip(_H_level_ns, _H_level_Ls, _H_level_binding_Es, _H_level_cross_sections)
end

#TODO broadcast over ν
function hydrogen_bound_absorption(ν, T, nH, nHe, ne, invU_H; 
                                   n_max=30, use_hubeny_generalization=false, taper=true)
    #TODO collapse some of these maps
    E_levels = map(1:n_max) do n
        RydbergH_eV - RydbergH_eV/n^2
    end
    ws = map(zip(1:n_max, E_levels)) do (n, E_level)
        hummer_mihalas_w(T, n, nH, nHe, ne; use_hubeny_generalization=use_hubeny_generalization)
    end

    partial_sum = 0.0
    for n = 1 : 4 # iterate over energy levels
        #the degeneracy is not included here
        ndens_state = ws[n] * exp(-E_levels[n] / (kboltz_eV * T))

        # TODO this is different from TS and Hillier
        # if the photon energy is greater than the ionization energy of the unperturbed atom with its
        # electron having primary quantum number n, the dissolved fraction is unity
        dissolved_fraction = if hplanck_eV * ν > RydbergH_eV/n^2
            1.0
        else
            # otherwise there is still some bf absorption because of level dissolution
            # this is the effective quantum number associated with the energy of the nth level plus 
            # that of the photon, i.e. the upper level
            n_eff = 1 / sqrt(1/n^2 - hplanck_eV*ν/RydbergH_eV)
            frac = 1 - hummer_mihalas_w(T, n_eff, nH, nHe, ne; use_hubeny_generalization=use_hubeny_generalization)
            if taper && n==1
                redcut = hplanck_eV * c_cgs / (RydbergH_eV * (1/n^2 - 1/(n+1)^2))
                λ = c_cgs / ν
                if λ > redcut
                    frac *= exp(-(c_cgs/ν - redcut)*1e6)
                end
            end
            frac
        end

        cross_section = 0
        for (n, L, Ebind, sigmas) in _H_cross_sections
            #Elow = sigmas.itp.knots[1][1]
            #println(round(hplanck_eV * c_cgs / Elow * 1e8))
            g = 2 * (2L + 1)
            cross_section += g * sigmas(hplanck_eV*ν)
        end

        partial_sum += ndens_state * cross_section * dissolved_fraction
    end
    #factor of 10^-18 converts cross-sections from megabarns to cm^2
    nH * invU_H * partial_sum * (1.0 - exp(-hplanck_eV * ν / (kboltz_eV * T))) * 1e-18
end

"""
    hummer_mihalas_w(T, n, E_level, nH, nHe, ne)


!!!note
    This is experimental, and not used by Korg for spectal synthesis.

Calculate the corection, w, to the occupation fraction of a hydrogen energy level using the 
occupation probability formalism from Hummer and Mihalas 1988, optionally with the generalization by 
Hubeny+ 1994.  (Sometimes Daepen+ 1987 is cited instead, but H&M seems to be where the theory 
originated. Presumably it was delayed in publication.)

The expression for w is in equation 4.71 of H&M.  K, the QM correction used in defined in equation 4.24.
Note that H&M's "N"s are numbers (not number densities), and their "V" is volume.  These quantities 
apear only in the form N/V, so we use the number densities instead.

This is based partially on Paul Barklem and Kjell Eriksson's 
[WCALC fortran routine](https://github.com/barklem/hlinop/blob/master/hbop.f) 
(part of HBOP.f), which is used by (at least) Turbospectrum and SME.  As in that routine, we do 
consider hydrogen and helium as the relevant neutral species, and assume them to be in the ground 
state.  All ions are assumed to have charge 1.  Unlike that routine, the generalization to the 
formalism from Hubeny+ 1994 is turned off by default because I haven't closely checked it.  The 
difference effects the charged_term only, and temerature is only used when 
`use_hubeny_generalization` is set to `true`.
"""
function hummer_mihalas_w(T, n_eff, nH, nHe, ne; use_hubeny_generalization=false)
    # contribution to w from neutral species (neutral H and He, in this implementation)
    # this is sqrt<r^2> assuming l=0.  I'm unclear why this is the approximation barklem uses.
    r_level = sqrt(5/2*n_eff^4 + 1/2*n_eff^2)*bohr_radius_cgs 
    # how do I reproduce this helium radius?
    neutral_term = nH * (r_level + sqrt(3)*bohr_radius_cgs)^3 + nHe * (r_level + 1.02bohr_radius_cgs)^3
        
    # contributions to w from ions (these are assumed to be all singly ionized, so n_ion = n_e)
    # K is a  QM correction defined in H&M '88 equation 4.24
    K = if n_eff > 3
        # WCALC drops the final factor, which is nearly within 1% of unity for all n
        16/3 * (n_eff/(n_eff+1))^2 * ((n_eff + 7/6)/(n_eff^2 + n_eff + 1/2))
    else
        1.0
    end
    χ = RydbergH_eV / n_eff^2 * eV_to_cgs # binding energy
    e = electron_charge_cgs
    charged_term = if use_hubeny_generalization
        # this is a straight line-by-line port from HBOP. Review and rewrite if used.
        if (ne > 10) && (T > 10) 
            A = 0.09 * exp(0.16667*log(ne)) / sqrt(T)
            X = exp(3.15*log(1+A))
            BETAC = 8.3e14 * exp(-0.66667*log(ne)) * K / n_eff^4
            F = 0.1402*X*BETAC^3 / (1+0.1285*X*BETAC*sqrt(BETAC))
            log(F/(1+F)) / (-4π/3)
        else
            0
        end
    else
        16 * ((e^2)/(χ * sqrt(K)))^3 * ne
    end
        
    exp(-4π/3 * (neutral_term + charged_term))
end


"""
    hummer_mihalas_U_H(T, nH, nHe, ne)

!!!note
    This is experimental, and not used by Korg for spectal synthesis.

Calculate the partition function of neutral hydrogen using the occupation probability formalism
from Hummer and Mihalas 1988.  See [`hummer_mihalas_w`](@ref) for details.
"""
function hummer_mihalas_U_H(T, nH, nHe, ne; use_hubeny_generalization=false)
    # Used by hummer_mihalas_U_H, below. These are from NIST, but it would be nice to generate them on the fly.
    hydrogen_energy_levels = [0.0, 10.19880615024, 10.19881052514816, 10.19885151459, 12.0874936591, 12.0874949611, 12.0875070783, 12.0875071004, 12.0875115582, 12.74853244632, 12.74853299663, 12.7485381084, 12.74853811674, 12.74853999753, 12.748539998, 12.7485409403, 13.054498182, 13.054498464, 13.054501074, 13.054501086, 13.054502042, 13.054502046336, 13.054502526, 13.054502529303, 13.054502819633, 13.22070146198, 13.22070162532, 13.22070313941, 13.22070314214, 13.220703699081, 13.22070369934, 13.220703978574, 13.220703979103, 13.220704146258, 13.220704146589, 13.220704258272, 13.320916647, 13.32091675, 13.320917703, 13.320917704, 13.320918056, 13.38596007869, 13.38596014765, 13.38596078636, 13.38596078751, 13.385961022639, 13.4305536, 13.430553648, 13.430554096, 13.430554098, 13.430554262, 13.462451058, 13.462451094, 13.46245141908, 13.462451421, 13.46245154007, 13.486051554, 13.486051581, 13.486051825, 13.486051827, 13.486051916, 13.504001658, 13.504001678, 13.50400186581, 13.504001867, 13.50400193582]
    hydrogen_energy_level_degeneracies = [2, 2, 2, 4, 2, 2, 4, 4, 6, 2, 2, 4, 4, 6, 6, 8, 2, 2, 4, 4, 6, 6, 8, 8, 10, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 2, 2, 4, 4, 6, 2, 2, 4, 4, 6, 2, 2, 4, 4, 6, 2, 2, 4, 4, 6, 2, 2, 4, 4, 6, 2, 2, 4, 4, 6]
    hydrogen_energy_level_n = [1, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12]

    # for each level calculate the corection, w, and add the term to U
    # the expression for w comes from Hummer and Mihalas 1988 equation 4.71 
    U = 0.0
    for (E, g, n) in zip(hydrogen_energy_levels, hydrogen_energy_level_degeneracies, hydrogen_energy_level_n)
        n_eff = sqrt(RydbergH_eV / (RydbergH_eV - E)) # times Z, which is 1 for hydrogen
        w = hummer_mihalas_w(T, n_eff, nH, nHe, ne; use_hubeny_generalization=use_hubeny_generalization)
        U += w * g*exp(-E / (kboltz_eV * T))
    end
    U
end