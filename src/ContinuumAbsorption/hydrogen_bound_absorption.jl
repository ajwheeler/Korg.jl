using Interpolations: LinearInterpolation, Line
using HDF5

# load hydrogen bf cross sections
const _H_cross_sections = let
    h5open(joinpath(_data_dir, "bf_cross-sections", 
                                         "individual_H_cross-sections.h5")) do f
        sigmas = map(eachcol(read(f["E"])), eachcol(read(f["sigma"]))) do Es, σs
            LinearInterpolation(Es, σs, extrapolation_bc=Line())
        end
        # use the cross sections for the first 6 energy levels only.
        # the binding energy for n=7 corresponds to ~45,000 Å 
        collect(zip(read(f["n"]), sigmas))
    end
end

"""


Because MHD level dissolution applied to the the Lyman series limit leads to inflated cross-sections
in the visible, we don't use MHD for bf absorption from n=1.  This can be overridden by setting
`use_MHD_for_Lyman=true`, in which case you will also want to set `taper=true`, which the same 
tapering of the cross-section as [HBOP](https://github.com/barklem/hlinop/blob/master/hbop.f) to fix 
the problem.

The `use_hubeny_generalization` keyword argument ennables the generalization of the MHD from 
Hubeny 1994. It is experimental and switched off by default.
"""
function H_I_bf(νs, T, nH, nHe, ne, invU_H; 
                                   n_upper_max=40, use_hubeny_generalization=false, taper=false, 
                                   use_MHD_for_Lyman=false)
    #TODO collapse some of these maps
    E_levels = map(1:n_upper_max) do n
        RydbergH_eV - RydbergH_eV/n^2
    end
    ws = map(1:n_upper_max) do n
        hummer_mihalas_w(T, n, nH, nHe, ne; use_hubeny_generalization=use_hubeny_generalization)
    end

    total_cross_section = zeros(length(νs))
    for (n, sigmas) in _H_cross_sections
        #the degeneracy is not included here
        ndens_state = ws[n] * exp(-E_levels[n] / (kboltz_eV * T))

        dissolved_fraction = map(νs) do ν
            if hplanck_eV * ν > RydbergH_eV/n^2
                # if the photon energy is greater than the ionization energy of the unperturbed atom 
                # with its electron having primary quantum number n, the dissolved fraction is one.
                1.0
            elseif !use_MHD_for_Lyman && n==1
                # don't use MHD for the Lyman series limit since it leads to inflated cross-sections
                # far red of the limit
                0.0 
            else  # account for bf absorption redward of the limit because of level dissolution
                # the effective quantum number associated with the energy of the nth level plus the 
                # photon energy
                n_eff = 1 / sqrt(1/n^2 - hplanck_eV*ν/RydbergH_eV)
                # this could probably be interpolated without much loss of accuracy
                w_upper = hummer_mihalas_w(T, n_eff, nH, nHe, ne; 
                                            use_hubeny_generalization=use_hubeny_generalization)
                # w_upper/w[n] is the prob that the upper level is dissolved given that the lower isn't
                frac = 1 - w_upper/ws[n] 
                if taper 
                    # taper of the cross-section past a certain  wavelength redward of the jump, as 
                    # is done in HBOP. (Not ennabled in Korg calls to this function.)
                    redcut = hplanck_eV * c_cgs / (RydbergH_eV * (1/n^2 - 1/(n+1)^2))
                    λ = c_cgs / ν
                    if λ > redcut
                        frac *= exp(-(λ - redcut)*1e6)
                    end
                end
                frac
            end
        end
        
        cross_section = sigmas.(hplanck_eV .* νs)
        
        total_cross_section .+= ndens_state .* cross_section .* dissolved_fraction
    end
    #factor of 10^-18 converts cross-sections from megabarns to cm^2
    @. nH * invU_H * total_cross_section * (1.0 - exp(-hplanck_eV * νs / (kboltz_eV * T))) * 1e-18
end