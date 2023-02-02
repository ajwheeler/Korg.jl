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
function get_inv_nK(mol, T, log_equilibrium_constants) 
    inv_nK = (kboltz_cgs*T)^(n_atoms(mol) - 1) / 10^log_equilibrium_constants[mol](log(T))
    # this weird construction handles the fact that autodiff will set the partial derivatives of 
    # inv_nK to NaNs in cases where inv_nK is 0.  In these cases (∂ inv_nK / ∂ whatever) should also 
    # be numerically 0, so dropping them does no harm.
    if inv_nK == 0
        0
    else
        inv_nK
    end
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
        inv_nK = get_inv_nK(mol, T, log_equilibrium_constants)
        els = get_atoms(mol.formula)
        number_densities[mol] = prod(number_densities[Species(Formula(el), 0)] for el in els) * inv_nK
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
    inv_nKs = get_inv_nK.(molecules, T, Ref(log_equilibrium_constants))
                                    
    #`residuals!` puts the residuals the system of molecular equilibrium equations in `F`
    #`x` is a vector containing the number density of the neutral species of each element
    function residuals!(F, x)
        # Don't allow negative number densities.  This is a trick to bound the possible values 
        # of x. Taking the log was less performant in tests.
        x = abs.(x) 
                                    
        # LHS: total number of atoms, RHS: first through third ionization states
        F .= atom_number_densities .- ion_factors .* x
        for (m, inv_nK) in zip(molecules, inv_nKs)
            els = get_atoms(m.formula)
            n_mol = prod(x[el] for el in els) * inv_nK
            # RHS: atoms which are part of molecules
            for el in els
                F[el] -= n_mol
            end
        end
    end
end