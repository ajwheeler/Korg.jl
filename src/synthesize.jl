using Interpolations: linear_interpolation
import .ContinuumAbsorption: total_continuum_absorption
using .RadiativeTransfer

"""
    synthesize(atm, linelist, A_X, (λ_start, λ_stop); kwargs... )
    synthesize(atm, linelist, A_X, wavelength_ranges; kwargs... )

Compute a synthetic spectrum.

# Arguments

  - `atm`: the model atmosphere (see [`read_model_atmosphere`](@ref))
  - `linelist`: A vector of [`Line`](@ref)s (see [`read_linelist`](@ref),
    [`get_APOGEE_DR17_linelist`](@ref), and [`get_VALD_solar_linelist`](@ref)).
  - `A_X`: a vector containing the A(X) abundances (log(X/H) + 12) for elements from hydrogen to
    uranium.  (see [`format_A_X`](@ref))
  - The wavelengths at which to synthesize the spectrum.  They can be specified either as a
    pair `(λstart, λstop)`, or as a list of pairs `[(λstart1, λstop1), (λstart2, λstop2), ...]`.
  - `λ_start`: the lower bound (in Å) of the region you wish to synthesize.
  - `λ_stop`: the upper bound (in Å) of the region you wish to synthesize.
  - `λ_step` (default: 0.01): the (approximate) step size to take (in Å).

# Returns

A named tuple with keys:

  - `flux`: the output spectrum
  - `cntm`: the continuum at each wavelength
  - `intensity`: the intensity at each wavelength and mu value, and possibly each layer in the model
    atmosphere, depending on the radiative transfer scheme.
  - `alpha`: the linear absorption coefficient at each wavelength and atmospheric layer a Matrix of
    size (layers x wavelengths)
  - `mu_grid`: a vector of tuples containing the μ values and weights used in the radiative transfer
    calculation. Can be controlled with the `mu_values` keyword argument.
  - `number_densities`: A dictionary mapping `Species` to vectors of number densities at each
    atmospheric layer
  - `electron_number_density`: the electron number density at each atmospheric layer
  - `wavelengths`: The vector of vacuum wavelengths (in Å) over which the synthesis was performed.
    If `air_wavelengths=true` this will not be the same as the input wavelengths.
  - `subspectra`: A vector of ranges which can be used to index into `flux` to extract the spectrum
    for each range provided in `wavelength_ranges`.  If you use the standard `λ_start`, `λ_stop`,
    `λ_step` arguments, this will be a vector containing only one range.

# Example

to synthesize a spectrum between 5000 Å and 5100 Å, with all metal abundances set to
0.5 dex less than the solar value except carbon, except carbon, which we set to [C/H]=-0.25:

```
atm = read_model_atmosphere("path/to/atmosphere.mod")
linelist = read_linelist("path/to/linelist.vald")
A_X = format_A_X(-0.5, Dict("C" => -0.25))
solution = synthesize(atm, linelist, A_X, 5000, 5100)
```

# Optional arguments:

  - `vmic` (default: 0) is the microturbulent velocity, ``\\xi``, in km/s.
  - `line_buffer` (default: 10): the farthest (in Å) any line can be from the provided wavelength range
    before it is discarded.  If the edge of your window is near a strong line, you may have to turn
    this up.
  - `cntm_step` (default 1): the distance (in Å) between point at which the continuum opacity is
    calculated.
  - `hydrogen_lines` (default: `true`): whether or not to include H lines in the synthesis.
  - `use_MHD_for_hydrogen_lines` (default: `true`): whether or not to use the MHD occupation
    probability formalism for hydrogen lines. (MHD is always used for hydrogen bound-free absorption.)
  - `hydrogen_line_window_size` (default: 150): the maximum distance (in Å) from each hydrogen line
    center at which to calculate its contribution to the total absorption coefficient.
  - `mu_values` (default: 20): the number of μ values at which to calculate the surface flux, or a
    vector of the specific values to use when doing transfer in spherical geometry. If `mu_points` is
    an integer, the values are chosen per Gauss-Legendre integration. If they are specified directly,
    the trapezoid rule is used for the astrophysical flux. The default values is sufficient for
    accuracy at the 10^-3 level. Note that if you are using the default radiative transfer scheme,
    with a plane-parallel model atmosphere, the integral over μ is exact, so this parameter has no
    effect. The points and weights are returned in the `mu_grid` field of the output.
  - `line_cutoff_threshold` (default: `3e-4`): the fraction of the continuum absorption coefficient
    at which line profiles are truncated.  This has major performance impacts, since line absorption
    calculations dominate more syntheses.  Turn it down for more precision at the expense of runtime.
    The default value should effect final spectra below the 10^-3 level.
  - `electron_number_density_warn_threshold` (default: `Inf`): if the relative difference between the
    calculated electron number density and the input electron number density is greater than this value,
    a warning is printed.  By default, this warning is suppress (threshold is `Inf`) because it is
    very easily raised in cases where it is of no observable consequence.
    See also `electron_number_density_warn_min_value`, below.
  - `electron_number_density_warn_min_value` (default: `1e-4`): The minimum value of the ratio of
    the electron number density to the total number density at which a warning is printed.
  - `return_cntm` (default: `true`): whether or not to return the continuum at each wavelength.  If
    this is false, `solution.cntm` will be `nothing`.
  - `ionization_energies`, a `Dict` mapping `Species` to their first three ionization energies,
    defaults to `Korg.ionization_energies`.
  - `partition_funcs`, a `Dict` mapping `Species` to partition functions (in terms of ln(T)). Defaults
    to data from Barklem & Collet 2016, `Korg.default_partition_funcs`.
  - `equilibrium_constants`, a `Dict` mapping `Species` representing diatomic molecules to the base-10
    log of their molecular equilibrium constants in partial pressure form.  Defaults to data from
    Barklem and Collet 2016, `Korg.default_log_equilibrium_constants`.
  - `use_chemical_equilibrium_from` (default: `nothing`): Takes another solution returned by
    `synthesize`. When provided, the chemical equilibrium solution will be taken from this object,
    rather than being recomputed. This is physically self-consistent only when the abundances, `A_X`,
    and model atmosphere, `atm`, are unchanged.
  - `molecular_cross_sections` (default: `[]`): A vector of precomputed molecular cross-sections. See
    [`MolecularCrossSection`](@ref) for how to generate these. If you are using the default radiative
    transfer scheme, your molecular cross-sections should cover 5000 Å only if your linelist does.
  - `tau_scheme` (default: "linear"): how to compute the optical depth.  Options are "linear" and
    "bezier" (testing only--not recommended).
  - `I_scheme` (default: `"linear_flux_only"`): how to compute the intensity.  Options are `"linear"`,
    `"linear_flux_only"`, and `"bezier"`.  `"linear_flux_only"` is the fastest, but does not return the
    intensity values anywhere except at the top of the atmosphere.  "linear" performs an equivalent
    calculation, but stores the intensity at every layer. `"bezier"` is for testing and not
    recommended.
  - use_CUDA TODO
  - `verbose` (default: `false`): Whether or not to print information about progress, etc.
"""
function synthesize(atm::ModelAtmosphere, linelist, A_X::AbstractVector{<:Real},
                    wavelength_params...;
                    vmic::Real=1.0, line_buffer::Real=10.0, cntm_step::Real=1.0,
                    air_wavelengths=false, hydrogen_lines=true, use_MHD_for_hydrogen_lines=true,
                    hydrogen_line_window_size=150, mu_values=20, line_cutoff_threshold=3e-4,
                    electron_number_density_warn_threshold=Inf,
                    electron_number_density_warn_min_value=1e-4, return_cntm=true,
                    I_scheme="linear_flux_only", tau_scheme="anchored",
                    ionization_energies=ionization_energies,
                    partition_funcs=default_partition_funcs,
                    log_equilibrium_constants=default_log_equilibrium_constants,
                    molecular_cross_sections=[], use_chemical_equilibrium_from=nothing,
                    use_CUDA=false, verbose=false)
    wls = Wavelengths(wavelength_params...; air_wavelengths=air_wavelengths)
    if air_wavelengths
        @warn "The air_wavelengths keyword argument is deprecated and will be removed in a future release. Korg.air_to_vacuum can be used to do the convertion, or you can create a Korg.Wavelengths with air_wavelengths=true and pass that to synthesize."
    end

    if use_MHD_for_hydrogen_lines && (wls[end] > 13_000 * 1e-8)
        @warn "if you are synthesizing at wavelengths longer than 15000 Å (e.g. for APOGEE), setting use_MHD_for_hydrogen_lines=false is recommended for the most accurate synthetic spectra. This behavior may become the default in Korg 1.0."
    end

    # Add wavelength bounds check (Rayleigh scattering limitation)
    # we should really have an upper bound as well
    min_allowed_wavelength = 1300.0 * 1e-8  # cm
    if first(wls) < min_allowed_wavelength
        # this restruction comes from the Rayleigh scattering calculations
        throw(ArgumentError("Requested wavelength range ($(wls)) " *
                            " extends blue-ward of 1300 Å, the lowerest allowed wavelength."))
    end

    # work in cm
    cntm_step *= 1e-8
    line_buffer *= 1e-8

    # wavelenths at which to calculate the continuum
    cntm_windows = map(eachwindow(wls)) do (λstart, λstop)
        (λstart - line_buffer - cntm_step, λstop + line_buffer + cntm_step)
    end
    cntm_windows, _ = merge_bounds(cntm_windows)
    cntm_wls = Wavelengths([w[1]:cntm_step:w[2] for w in cntm_windows])

    #sort the lines if necessary
    if !issorted(linelist; by=l -> l.wl)
        @warn "Linelist isn't sorted.  Sorting it, which may cause a significant delay."
        linelist = sort(linelist; by=l -> l.wl)
    end
    # sort linelist and remove lines far from the synthesis region
    # first just the ones needed for α5 (fall back to default if they aren't provided)
    if tau_scheme == "anchored"
        linelist5 = get_alpha_5000_linelist(linelist)
    end
    # now the ones for the synthesis
    linelist = filter_linelist(linelist, wls, line_buffer)

    if length(A_X) != MAX_ATOMIC_NUMBER || (A_X[1] != 12)
        throw(ArgumentError("A(H) must be a 92-element vector with A[1] == 12."))
    end

    abs_abundances = @. 10^(A_X - 12) # n(X) / n_tot
    abs_abundances ./= sum(abs_abundances) #normalize so that sum(n(X)/n_tot) = 1

    #float-like type general to handle dual numbers
    α_type = promote_type(eltype(atm.layers).parameters..., eltype(linelist).parameters...,
                          eltype(wls), typeof(vmic), typeof.(abs_abundances)...)
    #the absorption coefficient, α, for each wavelength and atmospheric layer
    α = Matrix{α_type}(undef, length(atm.layers), length(wls))
    # each layer's absorption at reference λ (5000 Å). This isn't used with the "anchored" τ scheme.
    α5 = Vector{α_type}(undef, length(atm.layers))
    triples = map(enumerate(atm.layers)) do (i, layer)
        nₑ, n_dict = if isnothing(use_chemical_equilibrium_from)
            chemical_equilibrium(layer.temp, layer.number_density,
                                 layer.electron_number_density,
                                 abs_abundances, ionization_energies,
                                 partition_funcs, log_equilibrium_constants;
                                 electron_number_density_warn_threshold=electron_number_density_warn_threshold,
                                 electron_number_density_warn_min_value=electron_number_density_warn_min_value)
        else
            let sol = use_chemical_equilibrium_from
                (sol.electron_number_density[i],
                 Dict(s => sol.number_densities[s][i] for s in keys(sol.number_densities)))
            end
        end

        α_cntm_vals = reverse(total_continuum_absorption(eachfreq(cntm_wls), layer.temp, nₑ, n_dict,
                                                         partition_funcs))
        α_cntm_layer = linear_interpolation(cntm_wls, α_cntm_vals)
        α[i, :] .= α_cntm_layer(wls)

        if tau_scheme == "anchored"
            α5[i] = total_continuum_absorption([c_cgs / 5e-5], layer.temp, nₑ, n_dict,
                                               partition_funcs)[1]
        end

        nₑ, n_dict, α_cntm_layer
    end
    nₑs = first.(triples)
    #put number densities in a dict of vectors, rather than a vector of dicts.
    n_dicts = getindex.(triples, 2)
    number_densities = Dict([spec => [n[spec] for n in n_dicts]
                             for spec in keys(n_dicts[1])
                             if spec != species"H III"])
    #vector of continuum-absorption interpolators
    α_cntm = last.(triples)

    # line contributions to α5
    if tau_scheme == "anchored"
        α_cntm_5 = [_ -> a for a in copy(α5)] # lambda per layer
        if use_CUDA
            line_absorption_cuda!(view(α5, :, 1), linelist5, Korg.Wavelengths([5000]),
                                  get_temps(atm), nₑs, number_densities, partition_funcs,
                                  vmic * 1e5, α_cntm_5; cutoff_threshold=line_cutoff_threshold)
        else
            line_absorption!(view(α5, :, 1), linelist5, Korg.Wavelengths([5000]), get_temps(atm),
                             nₑs, number_densities, partition_funcs, vmic * 1e5, α_cntm_5;
                             cutoff_threshold=line_cutoff_threshold)
        end
        interpolate_molecular_cross_sections!(view(α5, :, 1), molecular_cross_sections,
                                              Korg.Wavelengths([5000]),
                                              get_temps(atm), vmic, number_densities)
    end

    source_fn = blackbody.((l -> l.temp).(atm.layers), wls')
    cntm = nothing
    if return_cntm
        cntm, _, _, _ = RadiativeTransfer.radiative_transfer(atm, α, source_fn, mu_values; α_ref=α5,
                                                             I_scheme=I_scheme, τ_scheme=tau_scheme)
    end

    if hydrogen_lines
        for (i, (layer, n_dict, nₑ)) in enumerate(zip(atm.layers, n_dicts, nₑs))
            nH_I = n_dict[species"H_I"]
            nHe_I = n_dict[species"He_I"]
            U_H_I = partition_funcs[species"H_I"](log(layer.temp))
            hydrogen_line_absorption!(view(α, i, :), wls, layer.temp, nₑ, nH_I, nHe_I,
                                      U_H_I, vmic * 1e5,
                                      hydrogen_line_window_size * 1e-8;
                                      use_MHD=use_MHD_for_hydrogen_lines)
        end
    end

    if use_CUDA
        line_absorption_cuda!(α, linelist, wls, get_temps(atm), nₑs, number_densities,
                              partition_funcs, vmic * 1e5, α_cntm;
                              cutoff_threshold=line_cutoff_threshold, verbose=verbose)
    else
        line_absorption!(α, linelist, wls, get_temps(atm), nₑs, number_densities, partition_funcs,
                         vmic * 1e5, α_cntm; cutoff_threshold=line_cutoff_threshold,
                         verbose=verbose)
    end
    interpolate_molecular_cross_sections!(α, molecular_cross_sections, wls, get_temps(atm), vmic,
                                          number_densities)

    flux, intensity, μ_grid, μ_weights = RadiativeTransfer.radiative_transfer(atm, α, source_fn,
                                                                              mu_values; α_ref=α5,
                                                                              I_scheme=I_scheme,
                                                                              τ_scheme=tau_scheme)

    (flux=flux, cntm=cntm, intensity=intensity, alpha=α, mu_grid=collect(zip(μ_grid, μ_weights)),
     number_densities=number_densities, electron_number_density=nₑs,
     wavelengths=wls .* 1e8, subspectra=subspectrum_indices(wls))
end

"""
    filter_linelist(linelist, wls, line_buffer)

Return a new linelist containing only lines within the provided wavelength ranges.
"""
function filter_linelist(linelist, wls, line_buffer; warn_empty=true)
    nlines_before = length(linelist)

    last_line_index = 0 # we need to keep track of this across iterations to avoid double-counting lines.
    sub_ranges = map(eachwindow(wls)) do (λstart, λstop)
        first_line_index = searchsortedfirst(linelist, (; wl=λstart - line_buffer); by=l -> l.wl)
        # ensure we don't double-count lines.
        first_line_index = max(first_line_index, last_line_index + 1)

        last_line_index = searchsortedlast(linelist, (; wl=λstop + line_buffer); by=l -> l.wl)

        first_line_index:last_line_index
    end
    linelist = vcat((linelist[r] for r in sub_ranges)...)

    if nlines_before != 0 && length(linelist) == 0 && warn_empty
        @warn "The provided linelist was not empty, but none of the lines were within the provided wavelength range."
    end
    linelist
end

"""
    get_alpha_5000_linelist(linelist)

Arguments:

  - `linelist`: the synthesis linelist, which will be used if it covers the 5000 Å region.

Return a linelist which can be used to calculate the absorption at 5000 Å, which is required for
the standard radiative transfer scheme.  If the provided linelist doesn't contain lines near 5000 Å,
use the built-in one. (see [`_load_alpha_5000_linelist`](@ref))
"""
function get_alpha_5000_linelist(linelist)
    # start by getting the lines in the provided linelist which effect the synthesis at 5000 Å
    # use a 21 Å line buffer, which 1 Å bigger than the coverage of the fallback linelist.
    linelist5 = filter_linelist(linelist, Korg.Wavelengths(5000, 5000), 21e-8; warn_empty=false)
    # if there aren't any, use the built-in one
    if length(linelist5) == 0
        _alpha_5000_default_linelist
        # if there are some, but they don't actually cross over 5000 Å, use the built-in one where they
        # aren't present
    elseif linelist5[1].wl > 5e-5
        ll = filter(_alpha_5000_default_linelist) do line
            line.wl < linelist5[1].wl
        end
        [ll; linelist5]
    elseif linelist5[end].wl < 5e-5
        ll = filter(_alpha_5000_default_linelist) do line
            line.wl > linelist5[end].wl
        end
        [linelist5; ll]
    else # if the built-in lines span 5000 Å, use them
        linelist5
    end
end

"""
    blackbody(T, λ)

The value of the Planck blackbody function for temperature `T` at wavelength `λ` [cm].
"""
function blackbody(T, λ)
    h = hplanck_cgs
    c = c_cgs
    k = kboltz_cgs

    2 * h * c^2 / λ^5 * 1 / (exp(h * c / λ / k / T) - 1)
end
