include("Peach1970.jl")

"""
    positive_ion_ff_absorption!(α_out::Vector{Real}, ν::Real, T::Real, number_densities::Dict, 
                                ne::Real, departure_coefficients=Peach1970.departure_coefficients)

Computes the linear absorption coefficient (in cm⁻¹) for all free-free interactions involving
positively charged atomic species. Uses the provided departure coefficients when they are available,
and the uncorrected hydrogenic approximation when they are not.

# Arguments
- `ν`: frequency in Hz
- `T`: temperature in K
- `number_densities` is a `Dict` mapping each `Species` to its number density.
- `ne`: the number density of free electrons.
- `departure_coefficients` (optional, defaults to 
   `ContinuumAbsorption.Peach1970.departure_coefficients`: 
   a dictionary mapping species to the departure coefficients for the ff process it participates in 
   (e.g. `species"C II"` maps to the C III ff departure coefficients--see note below).  Departure
   coefficients should be callables taking temperature and (photon energy / RydbergH / Zeff^2).  

!!! note
    A free-free interaction is named as though the species interacting with the free electron had 
    one more bound electron (in other words it's named as though the free-electron and ion were 
    bound together). To make matters more confusing, some sources present the free-free linear 
    absorption coefficient's formula in terms of the number density of the species in the 
    interaction's name (by including the Saha equation in the formula).

    To compute the absorption for a given free-free interaction, this function accesses elements 
    from the dictionary arguments that are associated with the species that participates in the 
    interaction. For example:
    - Si I ff absorption: uses the number density of Si II (`number_densities[species"Si II"]`) and
      checks `departure_coefficients[species"Si II"]` for departure coefficients.
    - Si II ff absorption: uses the number density of Si III (`number_densities[species"Si III"]`)
      and checks `departure_coefficients[species"Si III"]` for departure coefficients.
"""
function positive_ion_ff_absorption!(α_out::AbstractVector{<:Real}, νs::AbstractVector{<:Real}, 
                                     T::Real, number_densities::Dict, ne::Real;
                                     departure_coefficients=Peach1970.departure_coefficients())
    if !(_gauntff_T_bounds.lower <= T <= _gauntff_T_bounds.upper)
        return #if T is outside the supported range return without changing α_out
    end
    #these are the freqs which are within the supported range
    idx = (c_cgs/_gauntff_λ_bounds.upper) .< νs .< (c_cgs/_gauntff_λ_bounds.lower)

    #under the hydrogenic approximation, all species with a given net charge, Z, share the same
    #absorption cross-section. To minimize the cost, accumulate the number densities of all species
    #(without departure coefficients) that have Z=1 or Z=2 before calling hydrogenic_ff_absorption
    ndens_Z1 = 0.0
    ndens_Z2 = 0.0

    for (spec,ndens) in number_densities
        if spec.charge <= 0
            # skip neutral species. They don't participate in ff interations.
            # While Korg doesn't track negatively charged ions as a separate species at the moment,
            # skip them too, in case that changes.
            continue 
        elseif spec in keys(departure_coefficients)
            D = departure_coefficients[spec]
            # photon energy in Rydberg*Zeff^2, see equation (5) in Peach 1967 (NOT Peach 1970)
            # https://articles.adsabs.harvard.edu/pdf/1967MmRAS..71....1P
            σs = @. νs[idx] / spec.charge^2 * (hplanck_eV / Rydberg_eV) 
            # add directly to α_out if there is a departure coefficient
            @. α_out[idx] += (hydrogenic_ff_absorption(νs[idx], T, spec.charge, 
                                                       number_densities[spec], ne)*(1+D(T, σs)))
        else
            #sum up contributions of hydrogenic ff coeffs, add them to α_out at the end
            if (spec.charge == 1)     # e.g. O II
                ndens_Z1 += ndens
            elseif (spec.charge == 2) # e.g. O III
                ndens_Z2 += ndens
            else
                error("triply+ ionized species not supported")
            end
        end
    end

    #add contributions from species for which we use the uncorrected hydrogenic approximation
    @. α_out[idx] += hydrogenic_ff_absorption(νs[idx], T, 1, ndens_Z1, ne)
    @. α_out[idx] += hydrogenic_ff_absorption(νs[idx], T, 2, ndens_Z2, ne)
    ;
end

