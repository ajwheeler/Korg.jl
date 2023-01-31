using NLsolve

"""
    saha_ion_weights(T, nₑ, atom, ionization_energies, partition_functions)

Returns `(wII, wIII)`, where `wII` is the ratio of singly ionized to neutral atoms of a given 
element, and `wIII` is the ration of doubly ionized to neutral atoms.

arguments:
- temperature `T` [K]
- electron number density `nₑ` [cm^-3]
- atom, the atomic number of the element 
- `ionization_energies` is a collection indexed by integers (e.g. a `Vector`) mappping elements' 
   atomic numbers to their first three ionization energies
- `partition_funcs` is a `Dict` mapping species to their partition functions
"""
function saha_ion_weights(T, nₑ, atom, ionization_energies, partition_funcs::Dict)
    χI, χII, χIII = ionization_energies[atom]
    atom = Formula(atom)
    UI = partition_funcs[Species(atom, 0)](log(T))
    UII = partition_funcs[Species(atom, 1)](log(T))

    k = kboltz_eV
    transU = translational_U(electron_mass_cgs, T)
    
    wII =  2.0/nₑ * (UII/UI) * transU * exp(-χI/(k*T))
    wIII = if atom == Formula(1) # hydrogen
        0.0
    else
        UIII = partition_funcs[Species(atom, 2)](log(T))
        (wII*2.0/nₑ * (UIII/UII) * transU * exp(-χII/(k*T)))
    end
    wII, wIII
end

"""
    translational_U(m, T)

The (inverse) contribution to the partition function from the free movement of a particle.
Used in the Saha equation.

arguments
- `m` is the particle mass
- `T` is the temperature in K
"""
function translational_U(m, T)
     k = kboltz_cgs
     h = hplanck_cgs
     (2π*m*k*T/h^2)^1.5
end

"""
    get_inv_nK(mol, T, log_equilibrium_constants)

Given a molecule, `mol`, a temperature, `T`, and a dictionary of log equilbrium constants in partial
pressure form, return the inverse equilibrium constant in number density form, i.e. `1/nK` where 
`nK = n(A)n(B)/n(AB)`.
"""
function get_log_nK(mol, T, log_equilibrium_constants) 
    log_equilibrium_constants[mol](log(T)) - (n_atoms(mol) - 1)*log10(kboltz_cgs*T) 
end

"""
    chemical_equilibrium(T, nₜ, nₑ, absolute_abundances, ionization_energies, 
                         partition_fns, log_equilibrium_constants; x0=nothing)

Iteratively solve for the number density of each species. Returns a `Dict` mapping species to number 
densities.

TODO
arguments:
- the temperature, `T`, in K
- the number density of non-electron particles `nₜ`
- the electron number density `nₑ`
- A Dict of `absolute_abundances`, N_X/N_total
- a Dict of ionization energies, `ionization_energies`.  The keys of act as a list of all atoms.
- a Dict of partition functions, `partition_fns`
- a Dict of log molecular equilibrium constants, `equilibrium_constants`, in partial pressure form. 
  The keys of `equilibrium_constants` act as a list of all molecules.

The system of equations is specified with the number densities of the neutral atoms as free 
parameters.  Each equation specifies the conservation of a particular species, e.g. (simplified)

    n(O) = n(CO) + n(OH) + n(O I) + n(O II) + n(O III).

In this equation:
- `n(O)`, the number density of oxygen atoms in any form comes `absolute_abundances` and the total
number density (supplied later)
- `n(O I)` is a free parameter.  The numerical solver is varying this to satisfy the system of 
  equations.
- `n(O II)`, and `n(O III)` come from the Saha (ionization) equation given `n(O I)`
- `n(CO)` and `n(OH)` come from the molecular equilibrium constants K, which are precomputed 
  over a range of temperatures. 

Equilibrium constants are defined in terms of partial pressures, so e.g.

    K(OH)  ==  (p(O) p(H)) / p(OH)  ==  (n(O) n(H)) / n(OH)) kT
"""
function chemical_equilibrium(T, nₜ, nₑ, absolute_abundances, ionization_energies, 
                              partition_fns, log_equilibrium_constants; x0=nothing)
    if x0 === nothing
        #compute good first guess by neglecting molecules
        x0 = map(1:MAX_ATOMIC_NUMBER) do atom
            wII, wIII =  saha_ion_weights(T, nₑ, atom, ionization_energies, partition_fns)
            nₜ*absolute_abundances[atom] / (1 + wII + wIII)
        end
    end

    #numerically solve for equlibrium.
    sol = nlsolve(chemical_equilibrium_equations(T, nₜ, nₑ, absolute_abundances, ionization_energies,
                                                 partition_fns, log_equilibrium_constants),
                  x0; iterations=1_000, store_trace=true, ftol=nₜ*1e-12, autodiff=:forward)
    if !sol.f_converged
        error("Molecular equlibrium unconverged. \n", sol)
    end

    # start with the neutral atomic species.  Only the absolute value of sol.zero is
    # necessarilly correct.
    number_densities = Dict(Species.(Formula.(1:MAX_ATOMIC_NUMBER), 0) .=> abs.(sol.zero))
    #now the ionized atomic species
    for a in 1:MAX_ATOMIC_NUMBER
        wII, wIII = saha_ion_weights(T, nₑ, a, ionization_energies, partition_fns)
        number_densities[Species(Formula(a), 1)] = wII  * number_densities[Species(Formula(a), 0)]
        number_densities[Species(Formula(a), 2)] = wIII * number_densities[Species(Formula(a), 0)]
    end
    #now the molecules
    for mol in keys(log_equilibrium_constants)
        log_nK = get_log_nK(mol, T, log_equilibrium_constants)
        element_log_ns = (log10(number_densities[Species(Formula(el), 0)]) for el in get_atoms(mol.formula))
        number_densities[mol] = 10^(sum(element_log_ns) - log_nK)
    end

    number_densities
end

function chemical_equilibrium_equations(T, nₜ, nₑ, absolute_abundances, ionization_energies, 
                                        partition_fns, log_equilibrium_constants)
    molecules = collect(keys(log_equilibrium_constants))
    atom_number_densities = nₜ .* absolute_abundances
                                    
    # ion_factors is a vector of ( n(X I) + n(X II)+ n(X III) ) / n(X I) for each element X
    ion_factors = map(1:MAX_ATOMIC_NUMBER) do elem
        wII, wIII = saha_ion_weights(T, nₑ, elem, ionization_energies, partition_fns)
        (1 + wII + wIII)
    end
    # precalculate equilibrium coefficients. Here, K is in terms of number density, not partial
    # pressure, unlike those in equilibrium_constants.
    log_nKs = get_log_nK.(molecules, T, Ref(log_equilibrium_constants))
                                    
    #`residuals!` puts the residuals the system of molecular equilibrium equations in `F`
    #`x` is a vector containing the number density of the neutral species of each element
    function residuals!(F, x)
        # Don't allow negative number densities.  This is a trick to bound the possible values 
        # of x. Taking the log was less performant in tests.
        x = abs.(x) 
                                    
        # LHS: total number of atoms, RHS: first through third ionization states
        F .= atom_number_densities .- ion_factors .* x
        for (m, log_nK) in zip(molecules, log_nKs)
            els = get_atoms(m.formula)
            n_mol = 10^(sum(log10(x[el]) for el in els) - log_nK)
            # RHS: atoms which are part of molecules
            for el in els
                F[el] -= n_mol
            end
        end
    end
end

###################################################################################
# The methods below are experimental, and not used by Korg for spectral synthesis #
###################################################################################

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
function hummer_mihalas_w(T, n, E_level, nH, nHe, ne; use_hubeny_generalization=false)
        # contribution to w from neutral species (neutral H and He, in this implementation)
        # effective quantum numer of the H energy level
        n_eff = sqrt(Korg.RydbergH_eV / (Korg.RydbergH_eV - E_level)) # times Z, which is 1 for hydrogen
        #this is sqrt<r^2> assuming l=0.  I'm unclear why this is the approximation barklem uses.
        r_level = sqrt(5/2*n_eff^4 + 1/2*n_eff^2)*Korg.bohr_radius_cgs 
        # how do I reproduce this helium radius?
        neutral_term = nH * (r_level + sqrt(3)*Korg.bohr_radius_cgs)^3 + nHe * (r_level + 1.02Korg.bohr_radius_cgs)^3
        
        # contributions to w from ions (these are assumed to be all singly ionized, so n_ion = n_e)
        # K is a  QM correction defined in H&M '88 equation 4.24
        K = if n > 3
            # WCALC drops the final factor, which is nearly within 1% of unity for all n
            16/3 * (n/(n+1))^2 * ((n + 7/6)/(n^2 + n + 1/2))
        else
            1.0
        end
        χ = (Korg.RydbergH_eV - E_level) * Korg.eV_to_cgs # binding energy
        e = Korg.electron_charge_cgs
        charged_term = if use_hubeny_generalization
            # this is a straight line-by-line port from HBOP. Review and rewrite if used.
            if (ne > 10) && (T > 10) 
                A = 0.09 * exp(0.16667*log(ne)) / sqrt(T)
                X = exp(3.15*log(1+A))
                BETAC = 8.3e14 * exp(-0.66667*log(ne)) * K / n^4
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
        w = hummer_mihalas_w(T, n, E, nH, nHe, ne; use_hubeny_generalization=use_hubeny_generalization)
        U += w * g*exp(-E / (Korg.kboltz_eV * T))
    end
    U
end