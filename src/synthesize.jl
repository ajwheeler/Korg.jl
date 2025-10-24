using Interpolations: linear_interpolation
import .ContinuumAbsorption: total_continuum_absorption
using .RadiativeTransfer

"""
    SynthesisResult

The result of a synthesis. Returned by [`synthesize`](@ref).

# Fields

  - `flux`: the output spectrum, in units of erg/s/cm^4/Å
  - `cntm`: the continuum at each wavelength, in the same units as `flux`. This is the same as the
    spectrum obtained with an empty linelist and with H lines turned off.
  - `intensity`: the intensity at each wavelength and mu value, and possibly each layer in the model
    atmosphere, depending on the radiative transfer scheme.
  - `alpha`: the linear absorption coefficient at each wavelength and atmospheric layer, a Matrix of
    size (layers × wavelengths)
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
"""
@kwdef struct SynthesisResult
    # specify container types to make debugging easier, but more precise typing would be better
    flux::Vector
    cntm::Union{Vector,Nothing}
    intensity::Array # can be either matrix or 3-tensor
    alpha::Matrix
    mu_grid::Vector{Tuple}
    number_densities::Dict{Species,Vector}
    electron_number_density::Vector
    wavelengths::Vector
    subspectra::Vector
end

"""
    synthesize(atm, linelist, A_X, (λ_start, λ_stop); kwargs... )

Compute a synthetic spectrum. Returns a [`SynthesisResult`](@ref).

# Arguments

  - `atm`: the model atmosphere (see [`interpolate_marcs`](@ref) and [`read_model_atmosphere`](@ref))
  - `linelist`: A vector of [`Line`](@ref)s (see [`read_linelist`](@ref),
    [`get_APOGEE_DR17_linelist`](@ref), [`get_GES_linelist`](@ref),
    [`get_GALAH_DR3_linelist`](@ref), and [`get_VALD_solar_linelist`](@ref)).
  - `A_X`: a vector containing the A(X) abundances (log(X/H) + 12) for elements from hydrogen to
    uranium.  [`format_A_X`](@ref) can be used to easily create this vector.
  - The wavelengths at which to synthesize the spectrum.  They can be specified either as a
    pair `(λstart, λstop)`, or as a list of pairs `[(λstart1, λstop1), (λstart2, λstop2), ...]`, or
    as any other valid arguments [described here](@ref wldocs)

# Example

To synthesize a spectrum between 5000 Å and 5100 Å, with all metal abundances set to
0.5 dex less than the solar value except carbon, which we set to [C/H]=-0.25:

```
linelist = Korg.read_linelist("path/to/linelist.vald")
A_X = format_A_X(-0.5, Dict("C" => -0.25))
atm = Korg.interpolate_marcs(5777, 4.44, A_X)
result = synthesize(atm, linelist, A_X, (5000, 5100))
```

# Optional arguments:

  - `vmic` (default: 1): the microturbulent velocity, ``\\xi``, in km/s. This can be either a scalar
    value or a vector of values, one for each atmospheric layer.
  - `line_buffer` (default: 10): the farthest distance (in Å) any line can be from the provided wavelength range
    before it is discarded. If the edge of your window is near a strong line, you may have to turn
    this up.
  - `cntm_step` (default: 1): the distance (in Å) between points at which the continuum opacity is
    calculated.
  - `hydrogen_lines` (default: `true`): whether or not to include H lines in the synthesis.
  - `use_MHD_for_hydrogen_lines`: whether or not to use the MHD occupation
    probability formalism for hydrogen lines. (MHD is always used for hydrogen bound-free absorption.)
    This is false by default when your last wavelength is > 13,000 Å, true otherwise (discussed in
    [Wheeler+ 2024](https://ui.adsabs.harvard.edu/abs/2023arXiv231019823W/abstract)).
  - `hydrogen_line_window_size` (default: 150): the maximum distance (in Å) from each hydrogen line
    center at which to calculate its contribution to the total absorption coefficient.
  - `mu_values` (default: 20): the number of μ values at which to calculate the surface flux, or a
    vector of the specific values to use when doing transfer in spherical geometry. If `mu_points` is
    an integer, the values are chosen per Gauss-Legendre integration. If they are specified directly,
    the trapezoid rule is used for the astrophysical flux. The default value is sufficient for
    accuracy at the 10^-3 level. Note that if you are using the default radiative transfer scheme
    with a plane-parallel model atmosphere, the integral over μ is exact, so this parameter has no
    effect. The points and weights are returned in the `mu_grid` field of the output.
  - `line_cutoff_threshold` (default: `3e-4`): the fraction of the continuum absorption coefficient
    at which line profiles are truncated. This has major performance impacts, since line absorption
    calculations dominate most syntheses. Turn it down for more precision at the expense of runtime.
    The default value should affect final spectra at by a factor of 10^-3 or less compared to no cutoff.
  - `electron_number_density_warn_threshold` (default: `Inf`): if the relative difference between the
    calculated electron number density and the input electron number density is greater than this value,
    a warning is printed. By default, this warning is suppressed (threshold is `Inf`) because it is
    very easily raised in cases where it is of no observable consequence.
    See also `electron_number_density_warn_min_value`, below.
  - `electron_number_density_warn_min_value` (default: `1e-4`): The minimum value of the ratio of
    the electron number density to the total number density at which a warning is printed.
  - `return_cntm` (default: `true`): whether or not to return the continuum at each wavelength.  If
    this is false, `solution.cntm` will be `nothing`.
  - `use_internal_reference_linelist` (default: `true`): whether or not to use the internal linelist
    for computing the opacity at the reference wavelength, which is used for radiative transfer when
    `tau_scheme` is "anchored". If this is false, `linelist` will be used instead if it overlaps with
    the reference wavelength (5000 Å for MARCS model atmospheres). This does not apply when Korg
    does not have a built-in reference linelist for the reference wavelength.
  - `ionization_energies`, a `Dict` mapping `Species` to their first three ionization energies,
    defaults to `Korg.ionization_energies`.
  - `partition_funcs`: a `Dict` mapping `Species` to partition functions (as functions of ln(T)).
    Defaults to data from Barklem & Collet 2016, `Korg.default_partition_funcs`.
  - `equilibrium_constants`: a `Dict` mapping `Species` representing diatomic molecules to the base-10
    log of their molecular equilibrium constants in partial pressure form. Defaults to data from
    Barklem and Collet 2016, `Korg.default_log_equilibrium_constants`.
  - `use_chemical_equilibrium_from` (default: `nothing`): Takes another solution returned by
    `synthesize`. When provided, the chemical equilibrium solution will be taken from this object,
    rather than being recomputed. This is physically self-consistent only when the abundances, `A_X`,
    and model atmosphere, `atm`, are unchanged.
  - `molecular_cross_sections` (default: `[]`): A vector of precomputed molecular cross-sections. See
    [`MolecularCrossSection`](@ref) for how to generate these. If you are using the default radiative
    transfer scheme and set `use_internal_reference_linelist=false`, your molecular cross-sections
    should cover the reference wavelength only if your linelist does.
  - `tau_scheme` (default: "linear"): how to compute the optical depth.  Options are "linear" and
    "bezier" (testing only--not recommended).
  - `I_scheme` (default: `"linear_flux_only"`): how to compute the intensity.  Options are `"linear"`,
    `"linear_flux_only"`, and `"bezier"`.  `"linear_flux_only"` is the fastest, but does not return the
    intensity values anywhere except at the top of the atmosphere.  "linear" performs an equivalent
    calculation, but stores the intensity at every layer. `"bezier"` is for testing and not
    recommended.
"""
function synthesize(atm::ModelAtmosphere, linelist, A_X::AbstractVector{<:Real},
                    wavelength_params...;
                    vmic=1.0,
                    line_buffer::Real=10.0,
                    cntm_step::Real=1.0,
                    hydrogen_lines=true,
                    use_MHD_for_hydrogen_lines::Union{Nothing,Bool}=nothing,
                    hydrogen_line_window_size=150,
                    mu_values=20,
                    line_cutoff_threshold=3e-4,
                    electron_number_density_warn_threshold=Inf,
                    electron_number_density_warn_min_value=1e-4,
                    return_cntm=true,
                    use_internal_reference_linelist=true,
                    I_scheme="linear_flux_only",
                    tau_scheme="anchored",
                    ionization_energies=ionization_energies,
                    partition_funcs=default_partition_funcs,
                    log_equilibrium_constants=default_log_equilibrium_constants,
                    molecular_cross_sections=[],
                    use_chemical_equilibrium_from=nothing,)::SynthesisResult
    wls = if length(wavelength_params) > 1
        @warn "Passing multiple wavelength parameters to `synthesize` is deprecated.  Package them in a tuple instead: synthesize(atm, linelist, A_X, (λ_start, λ_stop))"
        Wavelengths(wavelength_params)
    else
        Wavelengths(wavelength_params[1])
    end

    if isnothing(use_MHD_for_hydrogen_lines)
        use_MHD_for_hydrogen_lines = wls[end] < 13_000 * 1e-8
    end

    # Add wavelength bounds check (Rayleigh scattering limitation)
    # we should really have an upper bound as well
    min_allowed_wavelength = 1300.0 * 1e-8  # cm
    if first(wls) < min_allowed_wavelength
        # this restriction comes from the Rayleigh scattering calculations
        throw(ArgumentError("Requested wavelength range ($(wls)) " *
                            " extends blue-ward of 1300 Å, the lowerest allowed wavelength."))
    end

    # work in cm
    cntm_step *= 1e-8
    line_buffer *= 1e-8

    # wavelengths at which to calculate the continuum
    cntm_windows = map(eachwindow(wls)) do (λstart, λstop)
        (λstart - line_buffer - cntm_step, λstop + line_buffer + cntm_step)
    end
    cntm_windows, _ = merge_bounds(cntm_windows)
    cntm_wls = Wavelengths([w[1]:cntm_step:w[2] for w in cntm_windows])

    #sort the lines if necessary
    if !issorted(linelist; by=l -> l.wl)
        @warn "Linelist isn't sorted. Sorting it, which may cause a significant delay."
        linelist = sort(linelist; by=l -> l.wl)
    end
    # sort linelist and remove lines far from the synthesis region
    # first just the ones needed for α5 (fall back to default if they aren't provided)
    if tau_scheme == "anchored"
        linelist5 = get_reference_wavelength_linelist(linelist, atm.reference_wavelength;
                                                      use_internal_reference_linelist)
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
                          eltype(wls), eltype(vmic), typeof.(abs_abundances)...)
    #the absorption coefficient, α, for each wavelength and atmospheric layer
    α = Matrix{α_type}(undef, length(atm.layers), length(wls))
    # each layer's absorption at reference λ. This isn't used with the "anchored" τ scheme.
    α_ref = Vector{α_type}(undef, length(atm.layers))
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
            α_ref[i] = total_continuum_absorption([c_cgs / atm.reference_wavelength], layer.temp,
                                                  nₑ, n_dict,
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
        α_cntm_ref = [_ -> a for a in copy(α_ref)] # lambda per layer
        line_absorption!(view(α_ref, :, 1),
                         linelist5,
                         Korg.Wavelengths([atm.reference_wavelength * 1e8]),
                         get_temps(atm),
                         nₑs,
                         number_densities,
                         partition_funcs,
                         vmic * 1e5,
                         α_cntm_ref;
                         cutoff_threshold=line_cutoff_threshold)
        interpolate_molecular_cross_sections!(view(α_ref, :, 1),
                                              molecular_cross_sections,
                                              Korg.Wavelengths([atm.reference_wavelength * 1e8]),
                                              get_temps(atm),
                                              vmic,
                                              number_densities)
    end

    source_fn = blackbody.((l -> l.temp).(atm.layers), wls')
    cntm = nothing
    if return_cntm
        cntm, _, _, _ = RadiativeTransfer.radiative_transfer(atm, α, source_fn, mu_values;
                                                             α_ref=α_ref,
                                                             I_scheme=I_scheme, τ_scheme=tau_scheme)
    end

    if hydrogen_lines
        for (i, (layer, n_dict, nₑ)) in enumerate(zip(atm.layers, n_dicts, nₑs))
            nH_I = n_dict[species"H_I"]
            nHe_I = n_dict[species"He_I"]
            U_H_I = partition_funcs[species"H_I"](log(layer.temp))
            ξ = (isa(vmic, Number) ? vmic : vmic[i]) * 1e5
            hydrogen_line_absorption!(view(α, i, :), wls, layer.temp, nₑ, nH_I, nHe_I,
                                      U_H_I, ξ, hydrogen_line_window_size * 1e-8;
                                      use_MHD=use_MHD_for_hydrogen_lines)
        end
    end

    line_absorption!(α, linelist, wls, get_temps(atm), nₑs, number_densities, partition_funcs,
                     vmic * 1e5, α_cntm; cutoff_threshold=line_cutoff_threshold)
    interpolate_molecular_cross_sections!(α, molecular_cross_sections, wls, get_temps(atm), vmic,
                                          number_densities)

    flux, intensity, μ_grid, μ_weights = RadiativeTransfer.radiative_transfer(atm, α, source_fn,
                                                                              mu_values;
                                                                              α_ref, I_scheme,
                                                                              τ_scheme=tau_scheme)
    # convert from erg/s/cm^5 to erg/s/cm^4/Å.
    flux .*= 1e-8
    cntm .*= 1e-8

    SynthesisResult(; flux, cntm, intensity, alpha=α,
                    mu_grid=collect(zip(μ_grid, μ_weights)), number_densities,
                    electron_number_density=nₑs, wavelengths=wls .* 1e8,
                    subspectra=subspectrum_indices(wls))
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
    get_reference_wavelength_linelist(linelist, reference_wavelength)

Arguments:

  - `linelist`: the user-specified linelist.
  - `reference_wavelength`: the reference wavelength of the model atmosphere.

Return a linelist which can be used to calculate the absorption at `reference_wavelength`, which is
required for the standard radiative transfer scheme.  If the provided linelist doesn't contain lines
near `reference_wavelength`, use a built-in one. (see [`_load_alpha_5000_linelist`](@ref))

At present the only built-in fallback is for 5000 Å, which is what MARCS uses. For all other
reference wavelengths, the user-provided linelist must contain lines near the reference wavelength.
"""
function get_reference_wavelength_linelist(linelist, reference_wavelength;
                                           use_internal_reference_linelist)
    if reference_wavelength == 5e-5 && use_internal_reference_linelist
        return _alpha_5000_default_linelist
    end

    # start by getting the lines in the provided linelist which effect the synthesis at 5000 Å
    # use a 21 Å line buffer, which 1 Å bigger than the coverage of the fallback linelist.
    filtered_linelist = filter_linelist(linelist, Korg.Wavelengths((5000, 5000)), 21e-8;
                                        warn_empty=false)

    # handle non-5000 Å reference wavelengths. There's no built-in fallback for these.
    if reference_wavelength != 5e-5
        if length(filtered_linelist) == 0
            throw(ArgumentError("The provided linelist does not contain any lines near the reference wavelength of $reference_wavelength Å. Korg has an internal reference linelist for 5000 Å (the MARCS default), but not for other values."))
        else
            filtered_linelist
        end
    end

    # if there aren't any, use the built-in one
    if length(filtered_linelist) == 0
        _alpha_5000_default_linelist
        # if there are some, but they don't actually cross over 5000 Å, use the built-in one where they
        # aren't present
    elseif filtered_linelist[1].wl > 5e-5
        ll = filter(_alpha_5000_default_linelist) do line
            line.wl < filtered_linelist[1].wl
        end
        [ll; filtered_linelist]
    elseif filtered_linelist[end].wl < 5e-5
        ll = filter(_alpha_5000_default_linelist) do line
            line.wl > filtered_linelist[end].wl
        end
        [filtered_linelist; ll]
    else # if the built-in lines span 5000 Å, use them
        filtered_linelist
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

    2 * h * c^2 / λ^5 * 1 / expm1(h * c / λ / k / T)
end
