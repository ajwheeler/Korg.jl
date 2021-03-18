using SpecialFunctions: expint

"""
    synthesize(atm, linelist, wls, [metallicity], [alpha]; abundances=Dict())

Solve the transfer equation to get the resultant astrophysical flux at each wavelength.
- `metallicity`, i.e. [metals/H] is log_10 solar relative
- `alpha`, i.e. [alpha/H] is log_10 solar relative
- `abundances` are A(X) format, i.e. A(x) = log_10(n_X/n_H), where n_X is the number density of X.
"""
function synthesize(atm, linelist, wls::AbstractVector{Flt}, metallicity::Flt=0.0, 
                    alpha::Flt=0.00; abundances=Dict()) where Flt <: AbstractFloat
    #use solar abundances scaled by metallicity for those not provided
    elements = Set(get_elem(l.species) for l in linelist)
    for elem in elements
        if !(elem in keys(abundances))
            #I'm accessing the module global solar_abundances here, because this is a user-facing 
            #function.  Hopefully that's a reasonable line to draw.
            #TODO correctly decrease H (and He?) with increasing metallicity.  Currently this is 
            #TODO implement alpha
            #slightly broken
            abundances[elem] = solar_abundances[elem] * 10.0^metallicity
        else
            abundances[elem] = 10.0^(abundances[elem]-12.0)
        end
    end

    #the linear opacity coefficient, α, for each wavelength and atmospheric layer
    α = Matrix{Flt}(undef, length(atm), length(wls))
    for (i, layer) in enumerate(atm)
        α[i, :] = line_opacity(linelist, wls, layer.temp, layer.electron_density, 
                               layer.number_density, abundances, atomic_masses, partition_funcs,
                               ionization_energies)
    end

    #the thickness of each atmospheric layer 
    #TODO fix edge case
    Δs = [0 ; diff((l->l.colmass).(atm)) ./ (l->l.density).(atm)[2:end]]
    
    #optical depth at each layer at each wavelenth
    τ = cumsum(α .* Δs, dims=1)

    #this is not just the standard solution to the transfer equation
    #the exponential integral function expint captures the integral over the disk of the star to get
    #the emergent astrophysical flux.  I was made aware of this form of the solution, by Edmonds+
    #1969 (https://ui.adsabs.harvard.edu/abs/1969JQSRT...9.1427E/abstract), which presents a form of
    #the equation.  You can verify it by substituting the variable of integration in the exponential
    #integal, t, with mu=1/t.
    source_fn = blackbody.((l->l.temp).(atm), wls')
    flux = map(zip(eachcol(τ), eachcol(source_fn))) do (τ_λ, S_λ)
        trapezoid_rule(τ_λ, S_λ .* expint.(τ_λ, 2))
    end
    #return the solution, along with quantities usefull for debugging.  I'm ambivalent RE whether 
    #the API should look like this long-term.
    (flux=flux, alpha=α, tau=τ, source_fn=source_fn)
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
