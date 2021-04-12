using NLsolve

"""
   setup_ionization_energies([filename])

Parses the table of ionization energies and returns it as a dictionary mapping elements to
their ionization energies, `[χ₁, χ₂, χ₃]`.
"""
function setup_ionization_energies(fname=joinpath(_data_dir, 
                                                  "BarklemCollet2016-ionization_energies.dat"))
    open(fname, "r") do f
        d = Dict{String, Vector{Float64}}()
        for line in eachline(f)
            if line[1] != '#'        
                toks = split(strip(line))
                #the first token is the atomiz number, which we ignore
                d[toks[2]] = parse.(Float64, toks[3:end])
            end
        end
        d
    end
end

"""
Returns `(wII, wIII)`, where `wII` is the ratio of singly ionized to neutral atoms of a given 
element, and `wIII` is the ration of doubly ionized to neutral atoms.

arguments:
- temperature `T` [K]
- electron number density `nₑ` [cm^-3]
- atom, the atomic symbol of the element 
- `ionization_energies` is a Dict mapping elements to their first three ionization energies
- `partition_funcs` is a Dict mapping species to their partition functions
"""
function saha_ion_weights(T, nₑ, atom::String, ionization_energies::Dict, partition_funcs::Dict)
    χI, χII, χIII = ionization_energies[atom]
    UI = partition_funcs[atom*"_I"](T)
    UII = partition_funcs[atom*"_II"](T)

    k = kboltz_eV
    transU = translational_U(electron_mass_cgs, T)
    
    wII =  (2.0/nₑ * (UII/UI) * transU * exp(-χI/k/T))
    wIII = if atom == "H"
        0.0
    else
        UIII = partition_funcs[atom*"_III"](T)
        (wII*2.0/nₑ * (UIII/UII) * transU * exp(-χII/k/T))
    end
    wII, wIII
end
"convieience method for easier testing"
function saha_ion_weights(T, nₑ, χs::Vector{<:AbstractFloat}, Us::Vector{Function})
    ionization_energies = Dict(["X" => χs])
    partition_funcs = Dict(["X_I" => Us[1], "X_II" => Us[2], "X_III" => Us[3]])
    saha_ion_weights(T, nₑ, "X", ionization_energies, partition_funcs)
end

"""
The contribution to the partition function from the free movement of a particle.
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
Returns a NamedTuple representing the system of equations specifying molecular equilibrium.

arguments:
- A Dict of `absolute_abundances`, N_X/N_total
- a Dict of ionization energies, `ionization_energies`.  The keys of `ionization_energies` act as a 
  list of all atoms.
- a Dict of partition functions, `partition_fns`
- a Dict of log moleculare equilibrium constants, `equilibrium_constants`, in partial pressure form. 
  The keys of `equilibrium_constants` act as a list of all molecules.

The system of equations is specified with the number densities of the neutral atoms as free 
parameters.  Each equation specifies the conservation of a particular species, e.g. (simplified)

n(O) = n(CO) + n(OH) + n(O I) + n(O II) + n(O III).

In this equation:
- n(O), the number density of oxygen atoms in any form comes `absolute_abundances` and the total
number density (supplied later)
- n(O I) is a free parameter.  The numerical solver is varying this to satisfy the system of 
  equations.
- n(O II), and n(O III) come from the Saha (ionization) equation given n(O I)
- n(CO) and n(OH) come from the molecular equilibrium constants K, which are precomputed 
  over a range of temperatures. 

K is defined in terms of partial pressures, so e.g.
K(OH)  ==  (p(O) p(H)) / p(OH)  ==  (n(O) + n(H)) / n(OH) kT
Alternatively, this could be computed via a Saha equation involving the molecular partition function
and dissolution energy.
"""
function molecular_equilibrium_equations(absolute_abundances, ionization_energies, partition_fns, 
                                         equilibrium_constants)
    atoms = collect(keys(ionization_energies))
    molecules = filter(keys(equilibrium_constants)) do m
        !('+' in m || '-' in m)
    end

    #construct a mapping from each atom to an integer these will be the indeces of each number 
    #density in the system equations, the xᵢs in the vector for of the residual equation, F(x) = 0.
    var_indices = Dict{String, Int}()
    for (i, a) in enumerate(atoms)
        var_indices[a] = i
    end

    #the residuals of the molecular equilibrium equations parametrized by T, electron number density
    #and number densities of each element [cm^-3]
    function system(T, nₜ, nₑ)
        #`residuals!` puts the residuals the system of molecular equilibrium equations in `F`
        #`x` is a vector containing the number density of the neutral species of each element
        function residuals!(F, x)
            for (i, a) in enumerate(atoms)
                #LHS: total number of atoms
                F[i] = nₜ * absolute_abundances[a]
                wII, wIII = saha_ion_weights(T, nₑ, a, ionization_energies, partition_fns)
                #RHS: first through third ionization states
                F[i] -= (1 + wII + wIII) * x[i]
            end
            for m in molecules
                el1, el2 = get_atoms(m)
                i1 = var_indices[el1]
                i2 = var_indices[el2]
                nₘ = x[i1] * x[i2] * kboltz_cgs * T / 10^equilibrium_constants[m](T)
                #RHS: subtract atoms which are part of mollecules
                F[i1] -= nₘ
                F[i2] -= nₘ
            end
        end
    end

    #passing atoms and molecules might seem a little weird architecturally, but it's partly in 
    #anticipation if automatically culling the list of species considered by what's in the line list
    #in the future
    (atoms=atoms, molecules=molecules, equations=system)
end

"""
Iteratively solve for the number density of each species. Returns a Dict mapping species to number 
densities.

arguments:
- the system of molecular equilibrium equations `MEQs` (the thing returned by 
`molecular_equilibrium_equations`)
- the temperature `T`, 
- the number density of non-electron particles `nₜ`
- the electron number density `nₑ`
- a Dict of N_X/N_total, `absolute_abundances`
- a Dict of ionization energies, `ionization_energies`
- a Dict of partition functions, `partition_fns`
- a Dict of log molecular equilibrium constants, `equilibrium_constants`, in partial pressure form.
"""
function molecular_equilibrium(MEQs, T, nₜ, nₑ, absolute_abundances, ionization_energies, 
                               partition_fns, equilibrium_constants; 
                               x0=[nₜ*absolute_abundances[a]* 0.8 for a in MEQs.atoms]) :: Dict
    #numerically solve for equlibrium.  This uses finite difference rather that autodiff bacause 
    #it's fast enough this way and requires fewer deps, but enabling autodiff is trivial
    sol = nlsolve(MEQs.equations(T, nₜ, nₑ), x0; iterations=20, store_trace=true, ftol=nₜ * 1e-12, 
                  autodiff=:forward)
    if !sol.f_converged
        error("Mollecular equlibrium unconverged", sol, "\n", sol.trace)
    end

    #start with the neutral atomic species
    number_densities = Dict((MEQs.atoms .* "_I") .=> sol.zero)
    #now the ionized atomic species
    for a in MEQs.atoms
        wII, wIII = saha_ion_weights(T, nₑ, a, ionization_energies, partition_fns)
        number_densities[a*"_II"] = wII * number_densities[a*"_I"]
        number_densities[a*"_III"] = wIII * number_densities[a*"_I"]
    end
    #now the molecules
    for m in MEQs.molecules 
        el1, el2 = get_atoms(m)
        n₁ = number_densities[el1*"_I"]
        n₂ = number_densities[el2*"_I"]
        K = equilibrium_constants[m]
        number_densities[m*"_I"] = n₁ * n₂ * kboltz_cgs * T / 10^K(T)
    end

    number_densities
end
