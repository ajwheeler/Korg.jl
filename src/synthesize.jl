import Interpolations #for Interpolations.line
using Interpolations: LinearInterpolation
import ..ContinuumOpacity

"""
    synthesize(atm, linelist, λs, [metallicity, [alpha]]; abundances=Dict())

Solve the transfer equation in the model atmosphere `atm` with the transitions in `linelist` at the 
wavelengths `λs` [Å] to get the resultant astrophysical flux at each wavelength.

optional arguments:
- `metallicity`, i.e. [metals/H] is log_10 solar relative
- `vmic` (default: 0) is the microturbulent velocity, ξ, in km/s.
- `abundances` are A(X) format, i.e. A(x) = log_10(n_X/n_H), where n_X is the number density of X.
- `line_buffer` (default: 10): the farthest (in Å) any line can be from the provide wavelenth range 
   before it is discarded.  If the edge of your window is near a strong line, you may have to turn 
   this up.
- `cntm_step`: the wavelength resolution with which continuum opacities are calculated.
- `ionization_energies`, a Dict containing the first three ionization energies of each element, 
   defaults to `Korg.ionization_energies`.
- `partition_funcs`, a Dict mapping species to partition functions. Defaults to data from 
   Barklem & Collet 2016, `Korg.partition_funcs`.
- `equilibrium_constants`, a Dict mapping diatomic molecules to theirmolecular equilbrium constants
  in partial pressure form.  Defaults to data from Barklem and Collet 2016, 
  `Korg.equilibrium_constants`.

Uses solar abundances scaled by `metallicity` and for those not provided.
"""
function synthesize(atm, linelist, λs; metallicity::Real=0.0, vmic::Real=1.0, abundances=Dict(), 
                    line_buffer::Real=10.0, cntm_step::Real=1.0, 
                    ionization_energies=ionization_energies, partition_funcs=partition_funcs,
                    equilibrium_constants=equilibrium_constants)
    #work in cm
    λs = λs * 1e-8
    cntm_step *= 1e-8
    line_window *= 1e-8
    cntmλs = (λs[1] - line_window - cntm_step) : cntm_step : (λs[end] + line_window + cntm_step)

    #sort the lines if necessary and check that λs is sorted
    issorted(linelist; by=l->l.wl) || sort!(linelist, by=l->l.wl)
    if !issorted(λs)
        throw(ArgumentError("λs must be sorted"))
    end

    nlines = length(linelist)
    linelist = filter(l-> λs[1] - line_buffer*1e-8 <= l.wl <= λs[end] + line_buffer*1e-8, linelist)

    #impotent = setdiff(Set(keys(abundances)), elements)
    #if length(impotent) > 0
    #    @warn "Abundanc(es) for $(impotent) is specified but not in line list."
    #end

    abundances = get_absolute_abundances(atomic_symbols, metallicity, abundances)
    MEQs = molecular_equilibrium_equations(abundances, ionization_energies, partition_funcs, 
                                           equilibrium_constants)

    #the absorption coefficient, α, for each wavelength and atmospheric layer
    α_type = typeof(promote(atm[1].temp, length(linelist) > 0 ? linelist[1].wl : 1.0, λs[1], 
                            metallicity, vmic, abundances["H"])[1])
    α = Matrix{α_type}(undef, length(atm), length(λs))
    for (i, layer) in enumerate(atm)
        number_densities = molecular_equilibrium(MEQs, layer.temp, layer.number_density,
                                                 layer.electron_density)

        #Calculate the continuum absorption over cntmλs, which is a sparser grid, then construct an
        #interpolator that can be used to approximate it over a fine grid.
        α_cntm = LinearInterpolation(cntmλs,
                                    total_continuum_opacity(c_cgs ./ cntmλs, layer.temp, 
                                                            layer.electron_density, layer.density, 
                                                            number_densities, partition_funcs
                                                           ) * layer.density)
        α[i, :] = α_cntm.(λs)

        α[i, :] .+= line_absorption(linelist, λs, layer.temp, layer.electron_density, 
                                    number_densities, partition_funcs, vmic*1e5; α_cntm=α_cntm)

    end

    #the thickness of each atmospheric layer 
    Δcolmass = diff((l->l.colmass).(atm))
    Δs = 0.5([0 ; Δcolmass] + [Δcolmass; Δcolmass[end]]) ./ (l->l.density).(atm)

    τ = cumsum(α .* Δs, dims=1) #optical depth at each layer at each wavelenth

    source_fn = blackbody.((l->l.temp).(atm), λs')

    #The exponential integral function, expint, captures the integral over the disk of the star to 
    #get the emergent astrophysical flux. I was made aware of this form of the solution, by
    #Edmonds+ 1969 (https://ui.adsabs.harvard.edu/abs/1969JQSRT...9.1427E/abstract).
    #You can verify it by substituting the variable of integration in the exponential integal, t,
    #with mu=1/t.
    flux = map(zip(eachcol(τ), eachcol(source_fn))) do (τ_λ, S_λ)
        trapezoid_rule(τ_λ, S_λ .* exponential_integral_2.(τ_λ))
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
                           * "\"abundances\" to implicitly set the amount of H"))
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
The total continuum opacity, κ, at many frequencies, ν.

- `νs` are frequencies in Hz
- `T` is temperature in K
- `nₑ` is the electron number density in cm^-3
- `ρ` is the density in g cm^-3 
- `number_densities` is a `Dict` mapping each species to its number density
- `partition_funcs` is a `Dict mapping each species to its partition function
"""
function total_continuum_opacity(νs::Vector{F}, T::F, nₑ::F, ρ::F, number_densities::Dict, 
                                 partition_funcs::Dict) where F <: Real
    κ = zeros(F, length(νs))

    #TODO check all arguments

    #Hydrogen continuum opacities
    nH_I = number_densities["H_I"]
    nH_I_div_U = nH_I / partition_funcs["H_I"](T)
    κ += ContinuumOpacity.H_I_bf.(nH_I_div_U, νs, ρ, T) 
    κ += ContinuumOpacity.H_I_ff.(number_densities["H_II"], nₑ, νs, ρ, T)
    κ += ContinuumOpacity.Hminus_bf.(nH_I_div_U, nₑ, νs, ρ, T)
    κ += ContinuumOpacity.Hminus_ff.(nH_I_div_U, nₑ, νs, ρ, T)
    #κ += ContinuumOpacity.H2plus_bf_and_ff.(nH_I_div_U, number_densities["H_II"], νs, ρ, T)
    
    #He continuum opacities
    κ += ContinuumOpacity.He_II_bf.(number_densities["He_II"]/partition_funcs["He_II"](T), νs, ρ, T)
    #κ += ContinuumOpacity.He_II_ff.(number_densities["He_III"], nₑ, νs, ρ, T)
    # ContinuumOpacity.Heminus_ff is only valid for λ ≥ 5063 Å
    #κ += ContinuumOpacity.Heminus_ff.(number_densities["He_I"]/partition_funcs["He_I"](T),
    #                                  nₑ, νs, ρ, T)
    
    #electron scattering
    κ .+= ContinuumOpacity.electron_scattering(nₑ, ρ)
    
    κ
end


"""
    blackbody(T, λ)

The value of the Planck blackbody function for temperature `T` at wavelength `λ`.
"""
function blackbody(T, λ)
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

"""
Approximate second order exponential integral, E_2(x).  This stiches together several series 
expansions to get an approximation which is accurate within 1% for all `x`.
"""
function exponential_integral_2(x) 
    if x < 1.1
        _expint_small(x)
    elseif x < 2.5
        _expint_2(x)
    elseif x < 3.5
        _expint_3(x)
    elseif x < 4.5
        _expint_4(x)
    elseif x < 5.5
        _expint_5(x)
    elseif x < 6.5
        _expint_6(x)
    elseif x < 7.5
        _expint_7(x)
    elseif x < 9
        _expint_8(x)
    else
        _expint_large(x)
    end
end

function _expint_small(x) 
    #euler mascheroni constant
    ℇ = 0.57721566490153286060651209008240243104215933593992
    1 + ((log(x) + ℇ - 1) + (-0.5 + (0.08333333333333333 + (-0.013888888888888888 + 
                                                            0.0020833333333333333*x)*x)*x)*x)*x
end
function _expint_large(x)
    invx = 1/x
    exp(-x) * (1 + (-2 + (6 + (-24 + 120*invx)*invx)*invx)*invx)*invx
end
function _expint_2(x)
    x -= 2
    0.037534261820486914 + (-0.04890051070806112 + (0.033833820809153176 + (-0.016916910404576574 + 
                                          (0.007048712668573576 -0.0026785108140579598*x)*x)*x)*x)*x
end
function _expint_3(x)
    x -= 3
    0.010641925085272673   + (-0.013048381094197039   + (0.008297844727977323   + 
            (-0.003687930990212144   + (0.0013061422257001345  - 0.0003995258572729822*x)*x)*x)*x)*x
end
function _expint_4(x)
    x -= 4
    0.0031982292493385146  + (-0.0037793524098489054  + (0.0022894548610917728  + 
            (-0.0009539395254549051  + (0.00031003034577284415 - 8.466213288412284e-5*x )*x)*x)*x)*x
end
function _expint_5(x)
    x -= 5
    0.000996469042708825   + (-0.0011482955912753257  + (0.0006737946999085467  +
            (-0.00026951787996341863 + (8.310134632205409e-5   - 2.1202073223788938e-5*x)*x)*x)*x)*x
end
function _expint_6(x)
    x -= 6
    0.0003182574636904001  + (-0.0003600824521626587  + (0.00020656268138886323 + 
            (-8.032993165122457e-5   + (2.390771775334065e-5   - 5.8334831318151185e-6*x)*x)*x)*x)*x
end
function _expint_7(x)
    x -= 7
    0.00010350984428214624 + (-0.00011548173161033826 + (6.513442611103688e-5   + 
            (-2.4813114708966427e-5  + (7.200234178941151e-6   - 1.7027366981408086e-6*x)*x)*x)*x)*x
end
function _expint_8(x)
    x -= 8
    3.413764515111217e-5   + (-3.76656228439249e-5    + (2.096641424390699e-5   + 
            (-7.862405341465122e-6   + (2.2386015208338193e-6  - 5.173353514609864e-7*x )*x)*x)*x)*x
end

