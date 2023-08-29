using NLsolve, ForwardDiff

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
        wII*2.0/nₑ * (UIII/UII) * transU * exp(-χII/(k*T))
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

struct ChemicalEquilibriumError <: Exception
    msg::String
end
Base.showerror(io::IO, e::ChemicalEquilibriumError) = 
    print(io, "Chemical equillibrium failed: ",  e.msg)

"""
    chemical_equilibrium(T, nₜ, nₑ, absolute_abundances, ionization_energies, 
                         partition_fns, log_equilibrium_constants; x0=nothing)

Iteratively solve for the number density of each species. Returns a pair containing the electron 
number density and `Dict` mapping species to number densities.

arguments:
- the temperature, `T`, in K
- the number density of non-electron particles `nₜ`
- the electron number density `nₑ`
- A Dict of `absolute_abundances`, N_X/N_total
- a Dict of ionization energies, `ionization_energies`.  The keys of act as a list of all atoms.
- a Dict of partition functions, `partition_fns`
- a Dict of log molecular equilibrium constants, `log_equilibrium_constants`, in partial pressure form. 
  The keys of `equilibrium_constants` act as a list of all molecules.

keyword arguments:
- `x0` (default: `nothing`) is an initial guess for the solution (in the format internal to 
  `chemical_equilibrium`). If not supplied, a good guess is computed by neglecting molecules.
- `electron_number_density_warn_threshold` (default: `0.25`) is the fractional difference between 
  the calculated electron number density and the model atmosphere electron number density at which 
  a warning is issued.

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
function chemical_equilibrium(temp, nₜ, model_atm_nₑ, absolute_abundances, ionization_energies, 
                              partition_fns, log_equilibrium_constants;  
                              electron_number_density_warn_threshold=0.25)
    #compute good first guess by neglecting molecules
    neutral_fraction_guess = map(1:MAX_ATOMIC_NUMBER) do Z
        wII, wIII =  saha_ion_weights(temp, model_atm_nₑ, Z, ionization_energies, partition_fns)
        1 / (1 + wII + wIII) 
    end

    nₑ, neutral_fractions = solve_chemical_equilibrium(temp, nₜ, absolute_abundances, 
                                                       neutral_fraction_guess, model_atm_nₑ,
                                                       ionization_energies, partition_fns, 
                                                       log_equilibrium_constants)

    if abs((nₑ - model_atm_nₑ) / model_atm_nₑ) > electron_number_density_warn_threshold
        @warn "Electron number density differs from model atmosphere by a factor greater than $electron_number_density_warn_threshold. (calculated nₑ = $nₑ, model atmosphere nₑ = $model_atm_nₑ)"
    end 

    # start with the neutral atomic species.
    number_densities = Dict(Species.(Formula.(1:MAX_ATOMIC_NUMBER), 0) 
                            .=> nₜ .* absolute_abundances .* neutral_fractions)
    #now the ionized atomic species
    for a in 1:MAX_ATOMIC_NUMBER
        wII, wIII = saha_ion_weights(temp, nₑ, a, ionization_energies, partition_fns)
        number_densities[Species(Formula(a), 1)] = wII  * number_densities[Species(Formula(a), 0)]
        number_densities[Species(Formula(a), 2)] = wIII * number_densities[Species(Formula(a), 0)]
    end
    #now the molecules
    for mol in keys(log_equilibrium_constants)
        log_nK = get_log_nK(mol, temp, log_equilibrium_constants)
        element_log_ns = (log10(number_densities[Species(Formula(el), 0)]) for el in get_atoms(mol.formula))
        number_densities[mol] = 10^(sum(element_log_ns) - log_nK)
    end

    nₑ, number_densities
end

function solve_chemical_equilibrium(temp, nₜ, absolute_abundances, neutral_fraction_guess, nₑ_guess,
                                    ionization_energies, partition_fns, log_equilibrium_constants)
    zero = _solve_chemical_equilibrium(temp, nₜ, absolute_abundances, neutral_fraction_guess, nₑ_guess,
                                    ionization_energies, partition_fns, log_equilibrium_constants)
    nₑ = abs(zero[end]) * nₜ * 1e-5
    neutral_fractions = abs.(zero[1:end-1])
    nₑ, neutral_fractions 
end

function _solve_chemical_equilibrium(temp, nₜ, absolute_abundances, neutral_fraction_guess, nₑ_guess,
                                    ionization_energies, partition_fns, log_equilibrium_constants)
    #numerically solve for equlibrium.
    residuals! = setup_chemical_equilibrium_residuals(temp, nₜ, absolute_abundances, ionization_energies, 
                                                partition_fns, log_equilibrium_constants)

    x0 = [neutral_fraction_guess; nₑ_guess / nₜ * 1e5]

    # this wacky maneuver ensures that x0 has the approprate dual number type for autodiff
    # if that is going on.  I'm sure there's a better way...
    x0 = x0 .* (absolute_abundances[1] / absolute_abundances[1])

    sol = nlsolve(residuals!, x0; method=:newton, iterations=1_000, store_trace=true, ftol=1e-8, autodiff=:forward)

    if !sol.f_converged
        throw(ChemicalEquilibriumError("unconverged"))
    elseif !all(isfinite, sol.zero)
        throw(ChemicalEquilibriumError("solution contains non-finite values"))
    end

    sol.zero
end

# handle the case where a derivative is being taken with respect to T and ntot, but not abundances
function _solve_chemical_equilibrium(temp::ForwardDiff.Dual{T, V1, P}, 
                                    nₜ::ForwardDiff.Dual{T, V2, P},
                                    absolute_abundances::Vector{F},  # not duals!
                                    neutral_fraction_guess::Vector{ForwardDiff.Dual{T, V3, P}},
                                    nₑ_guess, # this type doesn't matter if the solver converges
                                    # Require that the types of the following be the same as the 
                                    # default, thus containing no duals. This is over-restrictive,
                                    # but I'm not sure how to enforce the more general condition
                                    # via the type system.
                                    ionization_energies::typeof(Korg.ionization_energies),
                                    partition_fns::typeof(Korg.default_partition_funcs),
                                    log_equilibrium_constants::typeof(Korg.default_log_equilibrium_constants)
                                    ) where {T, V1, V2, V3, P, F <: AbstractFloat}
    vtemp = ForwardDiff.value(temp)
    vnₜ = ForwardDiff.value(nₜ)
    vneutral_fraction_guess = ForwardDiff.value.(neutral_fraction_guess)
    vnₑ_guess = ForwardDiff.value(nₑ_guess)

    ptemp = ForwardDiff.partials(temp)
    pnₜ = ForwardDiff.partials(nₜ)
    partials = [ptemp pnₜ]'

    zero = _solve_chemical_equilibrium(vtemp, vnₜ, absolute_abundances, vneutral_fraction_guess, vnₑ_guess,
                                      ionization_energies, partition_fns, log_equilibrium_constants)

    residuals! = setup_chemical_equilibrium_residuals(vtemp, vnₜ, absolute_abundances, ionization_energies, 
                                                partition_fns, log_equilibrium_constants)

    tmp = similar(zero) # for storing results of residuals!. jacobian handles this nicely.
    drdx = ForwardDiff.jacobian((tmp, x) -> residuals!(tmp, x), tmp, zero)
    drdp = ForwardDiff.jacobian(tmp, [vtemp, vnₜ]) do tmp, p
        r! = setup_chemical_equilibrium_residuals(p[1], p[2], absolute_abundances, ionization_energies, 
                                            partition_fns, log_equilibrium_constants)
        r!(tmp, zero)
    end
    dxdp = -(drdx \ drdp)
    partial_zero = dxdp * partials

    dual_zero = map(zero, eachrow(partial_zero)) do v, p
        ForwardDiff.Dual{T}(v, p...)
    end

    dual_zero
end

function setup_chemical_equilibrium_residuals(T, nₜ, absolute_abundances, ionization_energies, 
                                        partition_fns, log_equilibrium_constants)
    molecules = collect(keys(log_equilibrium_constants))
    atom_number_densities = nₜ .* absolute_abundances
                                    
    # precalculate equilibrium coefficients. Here, K is in terms of number density, not partial
    # pressure, unlike those in equilibrium_constants.
    log_nKs = get_log_nK.(molecules, T, Ref(log_equilibrium_constants))

    # precompute the ratio of singly and doubly ionized to neutral atoms with factors of nₑ^-1 and 
    # nₑ^-2 divided out
    pairs = map(1:MAX_ATOMIC_NUMBER) do Z
        saha_ion_weights(T, 1, Z, ionization_energies, partition_fns)
    end
    wII_ne, wIII_ne2 = first.(pairs), last.(pairs)

    let nₜ=nₜ, log_nKs=log_nKs, molecules=molecules, atom_number_densities=atom_number_densities, wII_ne=wII_ne, wIII_ne2=wIII_ne2
        #`residuals!` puts the residuals the system of molecular equilibrium equations in `F`
        #`x` is a vector containing the number density of the neutral species of each element
        function residuals!(F, x)
            # the last is the electron number density in units of n_tot/10^5. The scaling allows us to 
            # better specify the tolerance for the solver.
            # Don't allow negative number densities.  This is a trick to bound the possible values 
            # of x. Taking the log was less performant in tests.
            nₑ = abs(x[end]) * nₜ * 1e-5 

            # the first 92 elements of x are the fraction of each element in it's neutral atomic form
            neutral_number_densities = atom_number_densities .* abs.(view(x, 1:MAX_ATOMIC_NUMBER))

            F[end] = 0

            # ion_factors is a vector of ( n(X I) + n(X II)+ n(X III) ) / n(X I) for each element X
            for Z in 1:MAX_ATOMIC_NUMBER
                wII = wII_ne[Z] / nₑ
                wIII = wIII_ne2[Z] / nₑ^2
                # LHS: total number of atoms, RHS: first through third ionization states
                F[Z] = atom_number_densities[Z] - (1 + wII + wIII) * neutral_number_densities[Z]
                # RHS: electrons freed from each ion
                F[end] += (wII + 2wIII) * neutral_number_densities[Z]
            end
            F[end] -= nₑ #LHS: total electron number density

            # now the first 92 elements of x are log10(neutral number densities)
            neutral_number_densities .= log10.(neutral_number_densities)
            for (m, log_nK) in zip(molecules, log_nKs)
                els = get_atoms(m.formula)
                n_mol = 10^(sum(neutral_number_densities[el] for el in els) - log_nK)
                # RHS: atoms which are part of molecules
                for el in els
                    F[el] -= n_mol
                end
            end

            F[1:end-1] ./= atom_number_densities
            F[end] /= nₑ * 1e-5
        end
    end
end


"""
    hummer_mihalas_w(T, n_eff, nH, nHe, ne; use_hubeny_generalization=false)

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

# hummer_mihalas_w is based partially on Paul Barklem and Nicolai Piskunov's HBOP routine. The 
# familly resemblance is limited to the "use_hubeny_generalization" option, which is not on be 
# defult, but we include license for HBOP here.
#
# Copyright (c) 2020, Paul Barklem and Nikolai Piskunov
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


"""
    hummer_mihalas_U_H(T, nH, nHe, ne)

!!!note
    This is experimental, and not used by Korg for spectral synthesis.

Calculate the partition function of neutral hydrogen using the occupation probability formalism
from Hummer and Mihalas 1988.  See [`hummer_mihalas_w`](@ref) for details.
"""
function hummer_mihalas_U_H(T, nH, nHe, ne; use_hubeny_generalization=false)
    # These are from NIST, but it would be nice to generate them on the fly.
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