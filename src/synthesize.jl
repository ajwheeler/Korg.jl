using Interpolations: LinearInterpolation
import ..ContinuumOpacity

"""
    synthesize(atm, linelist, λs; metallicity=0, abundances=Dict(), vmic=0, ... )

Solve the transfer equation in the model atmosphere `atm` with the transitions in `linelist` at the 
wavelengths `λs` [Å] to get the resultant astrophysical flux at each wavelength.

For efficiency reasons, `λs` must be an `AbstractRange`, such as `6000:0.01:6500`.  It can't be an 
arbitrary list of wavelengths.

Optional arguments:
- `metallicity`, i.e. [metals/H] is log_10 solar relative
- `abundances` is a `Dict` mapping atomic symbols to ``A(X)`` format abundances, i.e. 
   ``A(x) = \\log_{10}(n_X/n_\\mathrm{H}) + 12``, where ``n_X`` is the number density of ``X``.
   These override `metallicity`.
- `vmic` (default: 0) is the microturbulent velocity, ``\\xi``, in km/s.
- `line_buffer` (default: 10): the farthest (in Å) any line can be from the provide wavelenth range 
   before it is discarded.  If the edge of your window is near a strong line, you may have to turn 
   this up.
- `cntm_step` (default 1): the distance (in Å) between point at which the continuum opacity is 
  calculated.
- `hydrogen_lines` (default: `true`): whether or not to include H lines in the synthesis.
- `mu_grid`: the range of (surface) μ values at which to calculate the surface flux when doing 
   transfer in spherical geometry (when `atm` is a `ShellAtmosphere`).
- `line_cutoff_threshold` (default: `1e-3`): the fraction of the continuum absorption coefficient 
   at which line profiles are truncated.  This has major performance impacts, since line absorption
   calculations dominate more syntheses.  Turn it down for more precision at the expense of runtime.
   The default value should effect final spectra below the 10^-3 level.
- `ionization_energies`, a `Dict` mapping `Species` to their first three ionization energies, 
   defaults to `Korg.ionization_energies`.
- `partition_funcs`, a `Dict` mapping `Species` to partition functions. Defaults to data from 
   Barklem & Collet 2016, `Korg.partition_funcs`.
- `equilibrium_constants`, a `Dict` mapping `Species` representing diatomic molecules to their 
   molecular equilbrium constants in partial pressure form.  Defaults to data from 
   Barklem and Collet 2016, `Korg.equilibrium_constants`.
"""
function synthesize(atm::ModelAtmosphere, linelist, λs::AbstractRange; metallicity::Real=0.0, 
                    vmic::Real=1.0, abundances=Dict(), line_buffer::Real=10.0, cntm_step::Real=1.0, 
                    hydrogen_lines=true, mu_grid=0.05:0.05:1, line_cutoff_threshold=1e-3,
                    ionization_energies=ionization_energies, 
                    partition_funcs=partition_funcs, equilibrium_constants=equilibrium_constants)
    #work in cm
    λs = λs * 1e-8
    cntm_step *= 1e-8
    line_buffer *= 1e-8
    cntmλs = (λs[1] - line_buffer - cntm_step) : cntm_step : (λs[end] + line_buffer + cntm_step)

    #sort the lines if necessary and check that λs is sorted
    issorted(linelist; by=l->l.wl) || sort!(linelist, by=l->l.wl)
    if !issorted(λs)
        throw(ArgumentError("λs must be sorted"))
    end

    linelist = filter(l-> λs[1] - line_buffer*1e-8 <= l.wl <= λs[end] + line_buffer*1e-8, linelist)

    abundances = get_absolute_abundances(metallicity, abundances)
    MEQs = molecular_equilibrium_equations(abundances, ionization_energies, partition_funcs, 
                                           equilibrium_constants)
    n_dicts = map(atm.layers) do layer
         molecular_equilibrium(MEQs, layer.temp, layer.number_density, layer.electron_number_density)
    end
    number_densities = Dict([spec=>[n[spec] for n in n_dicts] for spec in keys(n_dicts[1])])

    #the absorption coefficient, α, for each wavelength and atmospheric layer
    α_type = typeof(promote(atm.layers[1].temp, length(linelist) > 0 ? linelist[1].wl : 1.0, λs[1], 
                            metallicity, vmic, abundances[1])[1])
    α = Matrix{α_type}(undef, length(atm.layers), length(λs))

    #Calculate the continuum absorption over cntmλs, which is a sparser grid, then construct an
    #interpolator that can be used to approximate it over a fine grid.
    α_cntm = map(zip(atm.layers, n_dicts)) do (layer, ns)
        LinearInterpolation(cntmλs, total_continuum_opacity(c_cgs ./ cntmλs, layer.temp, 
                                                            layer.electron_number_density, 
                                                            layer.density, ns, partition_funcs
                                                           ) * layer.density)
    end

    for (i, layer) in enumerate(atm.layers)
       if hydrogen_lines
           α[i, :] .= α_cntm[i].(λs)
           α[i, :] .+= hydrogen_line_absorption(λs, layer.temp, layer.electron_number_density, 
                                                number_densities[literals.H_I], 
                                                partition_funcs[literals.H_I], 
                                                hline_stark_profiles, vmic*1e5)
       end
    end

    line_absorption!(α, linelist, λs, [layer.temp for layer in atm.layers], 
                     [layer.electron_number_density for layer in atm.layers], number_densities,
                     partition_funcs, vmic*1e5; α_cntm=α_cntm, 
                     cutoff_threshold=line_cutoff_threshold)
    
    source_fn = blackbody.((l->l.temp).(atm.layers), λs')

    flux = if atm isa ShellAtmosphere
        rs = (l->l.r).(atm.layers)
        R = rs[1] + 0.5(rs[1] - rs[2])
        I = spherical_transfer(R, rs, α, source_fn, mu_grid) #μ-resolve intensity
        2π * [Korg.trapezoid_rule(mu_grid, mu_grid .* I) for I in eachrow(I)]
    else #atm isa PlanarAtmosphere
        #the thickness of each atmospheric layer 
        Δcolmass = diff((l->l.colmass).(atm.layers))
        Δs = 0.5([Δcolmass[1] ; Δcolmass] + [Δcolmass; Δcolmass[end]]) ./ (l->l.density).(atm.layers)
        τ = cumsum(α .* Δs, dims=1) #optical depth at each layer at each wavelenth

        map(zip(eachcol(τ), eachcol(source_fn))) do (τ_λ, S_λ)
            2π * transfer_integral(τ_λ, S_λ; plane_parallel=true)
        end
    end

    #idk whether we should return this extra stuff long-term, but it's useful for debugging
    (flux=flux, alpha=α, number_densities=number_densities)
end

"""
    get_absolute_abundances(metallicity, A_X)

Calculate ``n_X/n_\\mathrm{total}`` for each element X given some specified abundances, ``A(X)``.  Use the 
metallicity [``X``/H] to calculate those remaining from the solar values (except He).
"""
function get_absolute_abundances(metallicity, A_X::Dict) :: Vector{Number}
    if "H" in keys(A_X)
        throw(ArgumentError("A(H) set, but A(H) = 12 by definition. Adjust \"metallicity\" and "
                           * "\"abundances\" to implicitly set the amount of H"))
    end

    #populate dictionary of absolute abundaces
    abundances = map(0x01:Natoms) do elem
        if elem == 0x01 #hydrogen
            1.0
        elseif atomic_symbols[elem] in keys(A_X)
            10^(A_X[atomic_symbols[elem]] - 12.0)
        else
            #I'm accessing the module global solar_abundances here, but it doesn't make sense to 
            #make this an optional argument because this behavior can be completely overridden by 
            #specifying all abundances explicitely.
            Δ = (elem == 0x02 #= helium =#) ? 0.0 : metallicity
            10^(solar_abundances[elem] + Δ - 12.0)
        end
    end
    #now normalize so that sum(N_x/N_total) = 1
    abundances ./= sum(abundances)
    abundances
end

"""
    total_continuum_opacity(νs, T, nₑ, ρ, number_densities, partition_funcs)

The total continuum opacity, κ, at many frequencies, ν.

- `νs` are frequencies in Hz
- `T` is temperature in K
- `nₑ` is the electron number density in cm^-3
- `ρ` is the density in g cm^-3 
- `number_densities` is a `Dict` mapping each species to its number density
- `partition_funcs` is a `Dict` mapping each species to its partition function (e.g. 
  `Korg.partition_funcs`)
"""
function total_continuum_opacity(νs::Vector{F}, T::F, nₑ::F, ρ::F, number_densities::Dict, 
                                 partition_funcs::Dict) where F <: Real
    κ = zeros(F, length(νs))

    #TODO check all arguments

    #Hydrogen continuum opacities
    nH_I = number_densities[literals.H_I]
    nH_I_div_U = nH_I / partition_funcs[literals.H_I](T)
    κ += ContinuumOpacity.H_I_bf.(nH_I_div_U, νs, ρ, T) 
    κ += ContinuumOpacity.H_I_ff.(number_densities[literals.H_II], nₑ, νs, ρ, T)
    κ += ContinuumOpacity.Hminus_bf.(nH_I_div_U, nₑ, νs, ρ, T)
    κ += ContinuumOpacity.Hminus_ff.(nH_I_div_U, nₑ, νs, ρ, T)
    κ += ContinuumOpacity.H2plus_bf_and_ff.(nH_I_div_U, number_densities[literals.H_II], νs, ρ, T)
    
    #He continuum opacities
    κ += ContinuumOpacity.He_II_bf.(number_densities[literals.H_II] / 
                                    partition_funcs[literals.H_II](T), νs, ρ, T)
    #κ += ContinuumOpacity.He_II_ff.(number_densities[literals.He_III], nₑ, νs, ρ, T)
    # ContinuumOpacity.Heminus_ff is only valid for λ ≥ 5063 Å
    #κ += ContinuumOpacity.Heminus_ff.(number_densities[literals.He_I] / 
    #        partition_funcs[literals.He_I](T), nₑ, νs, ρ, T)
    
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
