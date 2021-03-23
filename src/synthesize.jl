using SpecialFunctions: expint
import ..ContinuumOpacity

"""
    synthesize(atm, linelist, λs, [metallicity], [alpha]; abundances=Dict())

Solve the transfer equation to get the resultant astrophysical flux at each wavelength.
- `metallicity`, i.e. [metals/H] is log_10 solar relative
- `abundances` are A(X) format, i.e. A(x) = log_10(n_X/n_H), where n_X is the number density of X.

Uses solar abundances scaled by `metallicity` and for those not provided.
"""
function synthesize(atm, linelist, λs::AbstractVector{F}, metallicity::F=0.0; abundances=Dict()
                   ) where F <: AbstractFloat
    #all the elements involved in either line or continuum opacity
    elements = Set(get_elem(l.species) for l in linelist) ∪ Set(["H", "He"])

    impotent = setdiff(Set(keys(abundances)), elements)
    if length(impotent) > 0
        @warn "Abundanc(es) for $(impotent) is specified but not in line list."
    end

    abundances = get_absolute_abundances(elements, metallicity, abundances)

    #TODO use `elements` below
    #the absorption coefficient, α, for each wavelength and atmospheric layer
    α = Matrix{F}(undef, length(atm), length(λs))
    nH = Vector{F}(undef, length(atm))
    for (i, layer) in enumerate(atm)
        number_densities = per_species_number_density(layer.number_density, layer.electron_density,
                                                      layer.temp, abundances)
        α[i, :] = line_opacity(linelist, λs, layer.temp, number_densities, atomic_masses, 
                               partition_funcs, ionization_energies)

        #continuum opacities are in frequency form, but this calculation is in wavelength form.
        νs = (c_cgs ./ λs) * 1e8
        #dλdν = (c_cgs ./ νs.^2) * 1e8
        α[i, :] += total_continuum_opacity(νs, layer.temp, layer.electron_density, layer.density,
                                           number_densities, partition_funcs) * layer.density
    end

    #the thickness of each atmospheric layer 
    #TODO fix edge case
    Δs = [0 ; diff((l->l.colmass).(atm)) ./ (l->l.density).(atm)[2:end]]
    τ = cumsum(α .* Δs, dims=1) #optical depth at each layer at each wavelenth

    source_fn = blackbody.((l->l.temp).(atm), λs')

    #This isn't the standard solution to the transfer equation.
    #The exponential integral function, expint, captures the integral over the disk of the star to 
    #get the emergent astrophysical flux. I was made aware of this form of the solution, by
    #Edmonds+ 1969 (https://ui.adsabs.harvard.edu/abs/1969JQSRT...9.1427E/abstract).
    #You can verify it by substituting the variable of integration in the exponential integal, t,
    #with mu=1/t.
    flux = map(zip(eachcol(τ), eachcol(source_fn))) do (τ_λ, S_λ)
        trapezoid_rule(τ_λ, S_λ .* expint.(2, τ_λ))
    end

    #return the solution, along with other quantities across wavelength and atmospheric layer.
    #idk whether we should return this extra stuff long-term, but it's useful for debugging
    (flux=flux, alpha=α, tau=τ, source_fn=source_fn)
end

"""
Calculate N_X/N_total for each X in `elements` given some `specified_abundances`, A(X).  Use the 
metallicity [X/H] to calculate those remaining from the solar values (except He).
"""
function get_absolute_abundances(elements, metallicity, A_X::Dict)::Dict
    if "H" in keys(A_X)
        throw(ArgumentError("A(H) set, but A(H) = 12 by definition. Adjust \"metallicity\" and "
                           * "\"abundnances\" to implicitely set the ammount of H"))
    end

    #populate dictionary of absolute abundaces
    abundances = Dict()
    for elem in elements
        if elem == "H"
            abundances[elem] = 1.0
        elseif elem in keys(A_X)
            abundances[elem] = 10^(A_X[elem] - 12.0)
        else
            #I'm accessing the module global solar_abundances here, but it doesn't make sense to 
            #make this an optional argument because this behavior can be completely overridden by 
            #specifying all abundances explicitely.
            Δ = elem == "He" ? 0.0 : metallicity
            abundances[elem] = 10^(solar_abundances[elem] + Δ - 12.0)
        end
    end
    #now normalize so that sum(N_x/N_total) = 1
    total = sum(values(abundances)) + sum([0; [10^(solar_abundances[elem] + metallicity - 12)
                                               for elem in atomic_symbols if ! (elem in elements)]])
    for elem in keys(abundances)
        abundances[elem] /= total
    end
    abundances
end

"""
Return a `Dict` which maps each species (i.e. Ba II) to number density [cm^-3]
- `nₜ` the number density of atoms and mollecules in cm^-3
- `nₑ` the number density of electrons in cm^-3
- `temp` the termperature in K
- `abundances` should be a `Dict` mapping species to absolute abundances between 0 and 1.
"""
function per_species_number_density(nₜ::Flt, nₑ::Flt, temp::Flt, abundances::Dict
                                   ) where Flt <: AbstractFloat
    #calculate the number density of each species
    n_densities = Dict()
    for (elem, abund) in abundances
        if elem == "H"
            #I'm not very happy with this, but I'm not sure what a better way to handle H is
            weights = saha(ionization_energies["H"], [partition_funcs["H_I"], 
                                                      partition_funcs["H_II"]],temp, nₑ)
            n_densities["H_I"] = nₜ * abund * weights[1]
            n_densities["H_II"] = nₜ * abund * weights[2]
        else
            weights = saha(ionization_energies[elem], 
                           [partition_funcs[elem * "_I"], partition_funcs[elem * "_II"], 
                              partition_funcs[elem * "_III"]], 
                           temp, nₑ)
            n_densities[elem * "_I"]   = nₜ * abund * weights[1]
            n_densities[elem * "_II"]  = nₜ * abund * weights[2]
            n_densities[elem * "_III"] = nₜ * abund * weights[3]
        end
    end
    n_densities
end

"""
The total continuum opacity, κ, at many frequencies, ν.

- `νs` are frequencies in Hz
- `T` is temperature in K
- `nₑ` is the electron number density in cm^-3
- `ρ` is the density in g cm^-3 
- `number_densities` is a `Dict` mapping each species to its number density
- `partition_funcs` is a `Dict mapping each species to its partition function
"""
function total_continuum_opacity(νs::Vector{F}, T::F, nₑ::F, ρ::F, number_densities::Dict, 
                                 partition_funcs::Dict) where F <: AbstractFloat
    κ = zeros(F, length(νs))

    #Hydrogen continuum opacities
    nH_I = number_densities["H_I"]
    nH_I_div_U = nH_I / partition_funcs["H_I"](T)
    κ += ContinuumOpacity.H_I_bf.(nH_I_div_U, νs, ρ, T) 
    κ += ContinuumOpacity.H_I_ff.(nH_I, nₑ, νs, ρ, T)
    κ += ContinuumOpacity.Hminus_bf.(nH_I_div_U, nₑ, νs, ρ, T)
    κ += ContinuumOpacity.Hminus_ff.(nH_I_div_U, nₑ, νs, ρ, T)
    #TODO uncomment when fixed
    #κ += ContinuumOpacity.H2plus_bf_and_ff.(nH_I_div_U, number_densities["H_II"], νs, ρ, T)
    
    #He continuum opacities
    nHe_II = number_densities["He_II"]
    nHe_II_div_U = nHe_II / partition_funcs["He_II"](T)
    κ += ContinuumOpacity.He_II_bf.(nHe_II_div_U, νs, ρ, T)
    #TODO fix
    #κ += ContinuumOpacity.He_II_ff.(nHe_II, nₑ, νs, ρ, T)
    κ += ContinuumOpacity.Heminus_ff.(number_densities["He_I"]/partition_funcs["He_I"](T), 
                                      nₑ, νs, ρ, T)
    
    #electron scattering
    κ .+= ContinuumOpacity.electron_scattering(nₑ, ρ)
    
    κ
end


"""
    blackbody(T, λ)

The value of the Planck blackbody function for temperature `T` at wavelength `λ`.
"""
function blackbody(T, λ)
    λ *= 1e-8 #convert to cm
    h = hplanck_cgs
    c = c_cgs
    k = kboltz_cgs

    2*h*c^2/λ^5 * 1/(exp(h*c/λ/k/T) - 1)
end

"""
    trapezoid_rule(xs, fs)

Approximate the integral from x₁ to x₂ of f(x) with the trapezoid rule given x-values `xs` and f(x)
values `fs`.

This should be good enough to numerically solve the transport equation, since model atmospheres
usually have carefully chosen knots.  We probably want to add higher-order aproximations later.
"""
function trapezoid_rule(xs, fs)
    Δs = diff(xs)
    weights = [0 ; Δs] + [Δs ; 0]
    sum(0.5 * weights .* fs)
end
