using Interpolations: LinearInterpolation
import ..ContinuumAbsorption: total_continuum_absorption

"""
    synthesize(atm, linelist, λ_start, λ_stop, [λ_step=0.01]; metallicity=0, abundances=Dict(), vmic=0, ... )

Solve the transfer equation in the model atmosphere `atm` with the transitions in `linelist` at the 
from `λ_start` to `λ_stop` in steps of `λ_step` (all Å) to get the resultant astrophysical flux at 
each wavelength.  An `AbstractRange` can also be probvided directly in place of `λ_start`, `λ_stop`,
and `λ_step`.

Returns a named tuple with keys:
- `flux`: the output spectrum
- `alpha`: the linear absorption coefficient at each wavelenth and atmospheric layer a Matrix of 
   size (layers x wavelengths)
- `number_densities`: A dictionary mapping `Species` to vectors of number densities at each 
   atmospheric layer
- `wavelengths`: The vacuum wavelenths (in Å) over which the synthesis was performed.  If 
  `air_wavelengths=true` this will not be the same as the input wavelenths.

## Optional arguments:
- `metallicity`, i.e. [metals/H] is the ``\\log_{10}`` solar-relative abundance of elements heavier 
   than He. It is overriden by `abundances`.
- `abundances` is a `Dict` mapping atomic symbols to [``X``/H] abundances.  (Set 
  `solar_relative=false` to use ``A(X)`` abundances instead.) These override `metallicity`.
- `vmic` (default: 0) is the microturbulent velocity, ``\\xi``, in km/s.
- `solar_relative` (default: true): When true, interpret abundances as being in \\[``X``/H\\] 
  (``\\log_{10}`` solar-relative) format.  When false, interpret them as ``A(X)`` abundances, i.e. 
   ``A(x) = \\log_{10}(n_X/n_\\mathrm{H}) + 12``, where ``n_X`` is the number density of ``X``.
   Note that abundances not specified default to the solar value, adjusted with `metallicity`, in
   either case.
- `solar_abundances` (default: `Korg.asplund_2020_solar_abundances`) is the set of solar abundances to 
  use, as a vector indexed by atomic number.  `Korg.asplund_2009_solar_abundances` and 
  `Korg.grevesse_2007_solar_abundances` are also provided for convienience.
- `air_wavelengths` (default: `false`): Whether or not the input wavelengths are air wavelenths to 
   be converted to vacuum wavelengths by Korg.  The conversion will not be exact, so that the 
   wavelenth range can internally be represented by an evenly-spaced range.  If the approximation 
   error is greater than `wavelength_conversion_warn_threshold`, an error will be thrown. (To do 
   wavelength conversions yourself, see [`air_to_vacuum`](@ref) and [`vacuum_to_air`](@ref).)
- `wavelength_conversion_warn_threshold` (default: 1e-4): see `air_wavelengths`.
- `line_buffer` (default: 10): the farthest (in Å) any line can be from the provided wavelenth range 
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
- `partition_funcs`, a `Dict` mapping `Species` to partition functions (in terms of ln(T)). Defaults 
   to data from Barklem & Collet 2016, `Korg.partition_funcs`.
- `equilibrium_constants`, a `Dict` mapping `Species` representing diatomic molecules to their 
   molecular equilbrium constants in partial pressure form.  Defaults to data from 
   Barklem and Collet 2016, `Korg.equilibrium_constants`.
"""
function synthesize(atm::ModelAtmosphere, linelist, λ_start, λ_stop, λ_step=0.01
                    ; air_wavelengths=false, wavelength_conversion_warn_threshold=1e-4, kwargs...)
    wls = if air_wavelengths
        len = Int(round((λ_stop - λ_start)/λ_step))+1
        vac_start, vac_stop = air_to_vacuum.((λ_start, λ_stop))
        vac_step = (vac_stop - vac_start) / (len-1)
        wls = StepRangeLen(vac_start, vac_step, len)
        max_diff = maximum(abs.(wls .- air_to_vacuum.(λ_start:λ_step:λ_stop)))
        if max_diff > wavelength_conversion_warn_threshold
            throw(ArgumentError("A linear air wavelength range can't be approximated exactly with a"
                                *"linear vacuum wavelength range. This solution differs by up to " * 
                                "$max_diff Å.  Adjust wavelength_conversion_warn_threshold if you" *
                                "want to suppress this error."))
        end
        wls
    else
        StepRangeLen(λ_start, λ_step, Int(round((λ_stop - λ_start)/λ_step))+1)
    end
    synthesize(atm, linelist, wls; kwargs...)
end
function synthesize(atm::ModelAtmosphere, linelist, λs::AbstractRange; metallicity::Real=0.0, 
                    vmic::Real=1.0, abundances=Dict(), line_buffer::Real=10.0, cntm_step::Real=1.0, 
                    hydrogen_lines=true, mu_grid=0:0.05:1, line_cutoff_threshold=1e-3,
                    solar_abundances=asplund_2020_solar_abundances, solar_relative=true,
                    ionization_energies=ionization_energies, 
                    partition_funcs=partition_funcs, equilibrium_constants=equilibrium_constants)
    #work in cm
    λs = λs * 1e-8
    cntm_step *= 1e-8
    line_buffer *= 1e-8
    cntmλs = (λs[1] - line_buffer - cntm_step) : cntm_step : (λs[end] + line_buffer + cntm_step)

    #sort the lines if necessary and check that λs is sorted
    issorted(linelist; by=l->l.wl) || sort!(linelist, by=l->l.wl)
    if step(λs) < 0
        throw(ArgumentError("λs must be in increasing order."))
    end

    #discard lines far from the wavelength range being synthesized
    linelist = filter(l-> λs[1] - line_buffer*1e-8 <= l.wl <= λs[end] + line_buffer*1e-8, linelist)

    abundances = get_absolute_abundances(metallicity, abundances, solar_abundances, solar_relative)
    MEQs = molecular_equilibrium_equations(abundances, ionization_energies, partition_funcs, 
                                           equilibrium_constants)

    sorted_cntmνs = c_cgs ./ reverse(cntmλs) #frequencies at which to calculate the continuum

    #float-like type general to handle dual numbers
    α_type = typeof(promote(atm.layers[1].temp, length(linelist) > 0 ? linelist[1].wl : 1.0, λs[1], 
                            metallicity, vmic, abundances[1])[1])
    #the absorption coefficient, α, for each wavelength and atmospheric layer
    α = Matrix{α_type}(undef, length(atm.layers), length(λs))
    α_cntm = Vector(undef, length(atm.layers))     #vector of continuum-absorption interpolators
    α5 = Vector{α_type}(undef, length(atm.layers)) #each layer's absorption at reference λ (5000 Å)
    n_dicts = Vector(undef, length(atm.layers))    #vector of (species -> number density) Dicts
    for (i, layer) in enumerate(atm.layers)
        n_dicts[i] = molecular_equilibrium(MEQs, layer.temp, layer.number_density, 
                                           layer.electron_number_density)

        α_cntm_vals = reverse(total_continuum_absorption(sorted_cntmνs, layer.temp,
                                                         layer.electron_number_density,
                                                         n_dicts[i], partition_funcs))
        α_cntm[i] = LinearInterpolation(cntmλs, α_cntm_vals)
        α[i, :] .= α_cntm[i].(λs)

        α5[i] = total_continuum_absorption([c_cgs/5e-5], layer.temp, layer.electron_number_density,
                                           n_dicts[i], partition_funcs)[1]

        if hydrogen_lines
            hydrogen_line_absorption!(view(α, i, :), λs, layer.temp, layer.electron_number_density, 
                                      n_dicts[i][species"H_I"], 
                                      partition_funcs[species"H_I"](log(layer.temp)), 
                                      hline_stark_profiles, vmic*1e5)
        end
    end

    #put number densities in a dict of vectors, rather than a vector of dicts.
    number_densities = Dict([spec=>[n[spec] for n in n_dicts] for spec in keys(n_dicts[1])])

    #add contribution of line absorption to α
    line_absorption!(α, linelist, λs, [layer.temp for layer in atm.layers], 
                     [layer.electron_number_density for layer in atm.layers], number_densities,
                     partition_funcs, vmic*1e5; α_cntm=α_cntm, 
                     cutoff_threshold=line_cutoff_threshold)

    source_fn = blackbody.((l->l.temp).(atm.layers), λs')
    flux = radiative_transfer(atm, α, source_fn, α5, mu_grid)

    (flux=flux, alpha=α, number_densities=number_densities, wavelengths=λs.*1e8)
end

"""
    get_absolute_abundances(metallicity, abundances, solar_abundances, solar_relative)

Calculate ``n_X/n_\\mathrm{total}`` for each element X given abundances in either [``X``/H] or 
``A(X)`` form.  See [`synthesize`](@ref) for a detailed description of the arguments.

arguments:
- `metallicity`: the 
"""
function get_absolute_abundances(metallicity::Real, abundances::Dict, 
                                 solar_abundances::AbstractVector, solar_relative::Bool
                                 )::Vector{Real}
    if "H" in keys(abundances)
        throw(ArgumentError("A(H) set, but A(H) = 12 by definition. Adjust \"metallicity\" and "
                           * "\"abundances\" to implicitly set the amount of H"))
    end

    #populate dictionary of absolute abundaces
    abundances = map(0x01:Natoms) do Z
        #abundances is indexed by atomic symbol, solar_abundances by atomic number
        elem = atomic_symbols[Z] 

        A_X = if elem in keys(abundances) #if explicitely set
            if solar_relative
                abundances[elem] + solar_abundances[Z]
            else
                abundances[elem]
            end
        else #if not set, use solar value adjusted for metallicity
            Δ = metallicity * (Z > 3) #only adjust for metals, not H or He
            solar_abundances[Z] + Δ
        end

        10^(A_X - 12)
    end
    
    abundances ./= sum(abundances) #normalize so that sum(N_x/N_total) = 1
    abundances
end

"""
    blackbody(T, λ)

The value of the Planck blackbody function for temperature `T` at wavelength `λ` [cm].
"""
function blackbody(T, λ)
    h = hplanck_cgs
    c = c_cgs
    k = kboltz_cgs

    2*h*c^2/λ^5 * 1/(exp(h*c/λ/k/T) - 1)
end
