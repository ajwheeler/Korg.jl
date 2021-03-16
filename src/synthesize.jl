using SpecialFunctions: expint

"""
    synthesize(TODO)

Solve the transfer equation to get the resultant astrophysical flux at each wavelength.
"""
function synthesize(atm, wavelength::Flt, metallicity::Flt=0.134, 
                    alpha=0.00; abundances=Dict()) where Flt <: AbstractFloat

    Δss = [0 ; diff((l->l.colmass).(atm)) ./ (l->l.density).(atm)[2:end]]

    αs = map(atm) do layer
        #calculate line and continuum opacities
        1e-8
    end
    τs = cumsum(αs .* Δss)
    source_fn = (l->blackbody(l.temp, wavelength)).(atm)

    #this is not just the standard solution to the transfer equation
    #the exponential integral function expint captures the integral over the disk of the star to get
    #the emergent astrophysical flux
    trapezoid_rule(τs, source_fn .* expint.(τs, 2))
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
