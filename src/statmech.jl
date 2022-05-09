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
    molecular_equilibrium_equations(absolute_abundances, ironization_energies, partiation_fns, equilibrium_constants)

Returns a NamedTuple representing the system of equations specifying molecular equilibrium.

arguments:
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
 
This could also be computed via a Saha equation involving the molecular partition function and 
dissolution energy.
"""
function molecular_equilibrium_equations(absolute_abundances, ionization_energies, partition_fns, 
                                         equilibrium_constants)
    atoms = 0x01:Natoms
    molecules = keys(equilibrium_constants)

    #the residuals of the molecular equilibrium equations parametrized by T, electron number density
    #and number densities of each element [cm^-3]
    function system(T, nₜ, nₑ)
        atom_number_densities = nₜ .* absolute_abundances
        ion_factors = map(atoms) do elem
            wII, wIII = saha_ion_weights(T, nₑ, elem, ionization_energies, partition_fns)
            (1 + wII + wIII)
        end
        #`residuals!` puts the residuals the system of molecular equilibrium equations in `F`
        #`x` is a vector containing the number density of the neutral species of each element
        function residuals!(F, x)
            #LHS: total number of atoms, RHS: first through third ionization states
            F .= atom_number_densities .- ion_factors .* x
            for m in molecules
                el1, el2 = get_atoms(m.formula)
                nₘ = x[el1] * x[el2] * kboltz_cgs * T / 10^equilibrium_constants[m](log(T))
                #RHS: subtract atoms which are part of mollecules
                F[el1] -= nₘ
                F[el2] -= nₘ
            end
        end
    end

    #passing atoms and molecules might seem a little weird architecturally, but it's partly in
    #anticipation of potentially culling the list of species considered 
    (atoms=atoms, molecules=molecules, equations=system, absolute_abundances=absolute_abundances,
     ionization_energies=ionization_energies, partition_fns=partition_fns, 
     equilibrium_constants=equilibrium_constants)
end

"""
    molecular_equilibrium(MEQS, T, nₜ, nₑ; x0)

Iteratively solve for the number density of each species. Returns a `Dict` mapping species to number 
densities.

arguments:
- the system of molecular equilibrium equations `MEQs` (the thing returned by 
`molecular_equilibrium_equations`)
- the temperature `T`, 
- the number density of non-electron particles `nₜ`
- the electron number density `nₑ`
- optionally, `x0`, a starting point for the solver
"""
function molecular_equilibrium(MEQs, T, nₜ, nₑ; x0=nothing)
    if x0 === nothing
        #compute good first guess by neglecting molecules
        x0 = map(MEQs.atoms) do atom
            wII, wIII =  saha_ion_weights(T, nₑ, atom, MEQs.ionization_energies, 
                                                  MEQs.partition_fns)
            nₜ*MEQs.absolute_abundances[atom] / (1 + wII + wIII)
        end
    end

    #numerically solve for equlibrium.
    sol = nlsolve(MEQs.equations(T, nₜ, nₑ), x0; iterations=200, store_trace=true, ftol=nₜ*1e-12,
                  autodiff=:forward)
    if !sol.f_converged
        error("Molecular equlibrium unconverged", sol, "\n", sol.trace)
    end

    #start with the neutral atomic species
    number_densities = Dict(Species.(Formula.(MEQs.atoms), 0) .=> sol.zero)
    #now the ionized atomic species
    for a in MEQs.atoms
        wII, wIII = saha_ion_weights(T, nₑ, a, MEQs.ionization_energies, MEQs.partition_fns)
        number_densities[Species(Formula(a), 1)] = wII  * number_densities[Species(Formula(a), 0)]
        number_densities[Species(Formula(a), 2)] = wIII * number_densities[Species(Formula(a), 0)]
    end
    #now the molecules
    for m in MEQs.molecules 
        el1, el2 = get_atoms(m.formula)
        n₁ = number_densities[Species(Formula(el1), 0)]
        n₂ = number_densities[Species(Formula(el2), 0)]
        logK = MEQs.equilibrium_constants[m]
        #add 1 to logK to conver to to cgs from mks
        number_densities[m] = n₁ * n₂ * kboltz_cgs * T / 10^logK(log(T))
    end

    number_densities
end
