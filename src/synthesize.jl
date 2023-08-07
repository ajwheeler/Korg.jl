using Interpolations: LinearInterpolation
import .ContinuumAbsorption: total_continuum_absorption
using .RadiativeTransfer

"""
    synthesize(atm, linelist, A_X, λ_start, λ_stop, [λ_step=0.01]; kwargs... )
    synthesize(atm, linelist, A_X, wavelength_ranges; kwargs... )

Compute a synthetic spectrum.

# Arguments
- `atm`: the model atmosphere (see [`read_model_atmosphere`](@ref))
- `linelist`: A vector of [`Line`](@ref)s (see [`read_linelist`](@ref), 
   [`get_APOGEE_DR17_linelist`](@ref), and [`get_VALD_solar_linelist`](@ref)).
- `A_X`: a vector containing the A(X) abundances (log(X/H) + 12) for elements from hydrogen to 
  uranium.  (see [`format_A_X`](@ref))
- `λ_start`: the lower bound (in Å) of the region you wish to synthesize.
- `λ_stop`: the upper bound (in Å) of the region you wish to synthesize.
- `λ_step` (default: 0.01): the (approximate) step size to take (in Å).

If you provide a vector of wavelength ranges in place of `λ_start` and `λ_stop`, the spectrum will 
be synthesized over each range with minimal overhead.
The ranges can be any Julia `AbstractRange`, for example: `5000:0.01:5010`.

# Returns 
A named tuple with keys:
- `flux`: the output spectrum
- `cntm`: the continuum at each wavelength
- `alpha`: the linear absorption coefficient at each wavelenth and atmospheric layer a Matrix of 
   size (layers x wavelengths)
- `number_densities`: A dictionary mapping species (as strings) to vectors of number densities at 
   each atmospheric layer
- `electron_number_density`: the electron number density at each atmospheric layer
- `wavelengths`: The vacuum wavelenths (in Å) over which the synthesis was performed.  If 
  `air_wavelengths=true` this will not be the same as the input wavelenths.
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
- `air_wavelengths` (default: `false`): Whether or not the input wavelengths are air wavelenths to 
   be converted to vacuum wavelengths by Korg.  The conversion will not be exact, so that the 
   wavelenth range can internally be represented by an evenly-spaced range.  If the approximation 
   error is greater than `wavelength_conversion_warn_threshold`, an error will be thrown. (To do 
   wavelength conversions yourself, see [`air_to_vacuum`](@ref) and [`vacuum_to_air`](@ref).)
- `wavelength_conversion_warn_threshold` (default: 1e-4): see `air_wavelengths`. (In Å.)
- `line_buffer` (default: 10): the farthest (in Å) any line can be from the provided wavelenth range 
   before it is discarded.  If the edge of your window is near a strong line, you may have to turn 
   this up.
- `cntm_step` (default 1): the distance (in Å) between point at which the continuum opacity is 
  calculated.
- `hydrogen_lines` (default: `true`): whether or not to include H lines in the synthesis.
- `use_MHD_for_hydrogen_lines` (default: `true`): whether or not to use the MHD occupation 
   probability formalism for hydrogen lines. (MHD is always used for hydrogen bound-free absorption.)
- `hydrogen_line_window_size` (default: 150): the mamximum distance (in Å) from each hydrogen line 
   center at which to calculate its contribution to the total absorption coefficient.
- `n_mu_points` (default: 20): the number of μ values at which to calculate the surface flux when doing 
   transfer in spherical geometry (when `atm` is a `ShellAtmosphere`). 20 points is sufficient for
   accuracy at the 10^-3 level.
- `line_cutoff_threshold` (default: `3e-4`): the fraction of the continuum absorption coefficient 
   at which line profiles are truncated.  This has major performance impacts, since line absorption
   calculations dominate more syntheses.  Turn it down for more precision at the expense of runtime.
   The default value should effect final spectra below the 10^-3 level.
- `electron_number_density_warn_threshold` (default: `1.0`): if the relative difference between the 
   calculated electron number density and the input electron number density is greater than this value,
   a warning is printed.  Set to `Inf` to suppress this warning.
- `return_cntm` (default: `true`): whether or not to return the continuum at each wavelength.  If 
   this is false, `solution.cntm` will be `nothing`.
- `ionization_energies`, a `Dict` mapping `Species` to their first three ionization energies, 
   defaults to `Korg.ionization_energies`.
- `partition_funcs`, a `Dict` mapping `Species` to partition functions (in terms of ln(T)). Defaults 
   to data from Barklem & Collet 2016, `Korg.default_partition_funcs`.
- `equilibrium_constants`, a `Dict` mapping `Species` representing diatomic molecules to the base-10
   log of their molecular equilbrium constants in partial pressure form.  Defaults to data from 
   Barklem and Collet 2016, `Korg.default_log_equilibrium_constants`.
- `bezier_radiative_transfer` (default: false): Use the radiative transfer scheme.  This is for 
   testing purposes only.
"""
function synthesize(atm::ModelAtmosphere, linelist, A_X, λ_start, λ_stop, λ_step=0.01; kwargs...)
    wls = [StepRangeLen(λ_start, λ_step, Int(round((λ_stop - λ_start)/λ_step))+1)]
    synthesize(atm, linelist, A_X, wls; kwargs...)
end
function synthesize(atm::ModelAtmosphere, linelist, A_X::AbstractVector{<:Real}, 
                    wl_ranges::AbstractVector{<:AbstractRange}; 
                    vmic::Real=1.0, line_buffer::Real=10.0, cntm_step::Real=1.0, 
                    air_wavelengths=false, wavelength_conversion_warn_threshold=1e-4,
                    hydrogen_lines=true, use_MHD_for_hydrogen_lines=true, 
                    hydrogen_line_window_size=150, n_mu_points=20, line_cutoff_threshold=3e-4, 
                    electron_number_density_warn_threshold=1.0, 
                    return_cntm=true,
                    bezier_radiative_transfer=false, ionization_energies=ionization_energies, 
                    partition_funcs=default_partition_funcs, 
                    log_equilibrium_constants=default_log_equilibrium_constants)

    # Convert air to vacuum wavelenths if necessary.
    if air_wavelengths
        wl_ranges = map(wl_ranges) do wls
            λ_start, λ_stop, λ_step = first(wls), last(wls), step(wls)
            len = Int(round((λ_stop - λ_start)/λ_step))+1
            vac_start, vac_stop = air_to_vacuum.((λ_start, λ_stop))
            vac_step = (vac_stop - vac_start) / (len-1)
            wls = StepRangeLen(vac_start, vac_step, len)
            max_diff = maximum(abs.(wls .- air_to_vacuum.(λ_start:λ_step:λ_stop)))
            if max_diff > wavelength_conversion_warn_threshold
                throw(ArgumentError("A linear air wavelength range can't be approximated exactly with a"
                                    *"linear vacuum wavelength range. This solution differs by up to " * 
                                    "$max_diff Å.  Adjust wavelength_conversion_warn_threshold if you "*
                                    "want to suppress this error."))
            end
            wls
        end
    end

    # work in cm
    wl_ranges = wl_ranges .* 1e-8 # broadbasting the = prevents wl_ranges from changing type
    cntm_step *= 1e-8
    line_buffer *= 1e-8

    # make sure wl_ranges are OK
    all_λs = vcat(wl_ranges...)
    if !issorted(all_λs) #TODO test
        throw(ArgumentError("wl_ranges must be sorted and non-overlapping"))
    end

    # wavelenths at which to calculate the continuum
    cntm_wl_ranges = map(wl_ranges) do λs 
        collect((λs[1] - line_buffer - cntm_step) : cntm_step : (λs[end] + line_buffer + cntm_step))
    end
    # eliminate portions where ranges overlap.  One fitting is merged, there will be functions for this.
    for i in 1:length(cntm_wl_ranges)-1 
        cntm_wl_ranges[i] = cntm_wl_ranges[i][cntm_wl_ranges[i] .< first(cntm_wl_ranges[i+1])]
    end
    cntmλs = vcat(cntm_wl_ranges...)
    # frequencies at which to calculate the continuum, as a single vector
    sorted_cntmνs = c_cgs ./ reverse(cntmλs) 

    #sort the lines if necessary
    issorted(linelist; by=l->l.wl) || sort!(linelist, by=l->l.wl)
    #discard lines far from the wavelength range being synthesized
    linelist = filter(linelist) do line
        map(wl_ranges) do wl_range
            wl_range[1] - line_buffer <= line.wl <= wl_range[end]
        end |> any
    end

    if length(A_X) != MAX_ATOMIC_NUMBER || (A_X[1] != 12)
        throw(ArgumentError("A(H) must be a 92-element vector with A[1] == 12."))
    end

    abs_abundances = @. 10^(A_X - 12) # n(X) / n_tot
    abs_abundances ./= sum(abs_abundances) #normalize so that sum(n(X)/n_tot) = 1

    #float-like type general to handle dual numbers
    α_type = promote_type(eltype(atm.layers).parameters..., eltype(linelist).parameters...,
                          eltype(all_λs), typeof(vmic), typeof.(abs_abundances)...)
    #the absorption coefficient, α, for each wavelength and atmospheric layer
    α = Matrix{α_type}(undef, length(atm.layers), length(all_λs))
    # each layer's absorption at reference λ (5000 Å)
    # This isn't used with bezier radiative transfer.
    α5 = Vector{α_type}(undef, length(atm.layers)) 
    triples = map(enumerate(atm.layers)) do (i, layer)
        nₑ, n_dict = chemical_equilibrium(layer.temp, layer.number_density, 
                                          layer.electron_number_density, 
                                          abs_abundances, ionization_energies, 
                                          partition_funcs, log_equilibrium_constants; 
                                          electron_number_density_warn_threshold=electron_number_density_warn_threshold)

        α_cntm_vals = reverse(total_continuum_absorption(sorted_cntmνs, layer.temp, nₑ, n_dict, 
                                                         partition_funcs))
        α_cntm_layer = LinearInterpolation(cntmλs, α_cntm_vals)
        α[i, :] .= α_cntm_layer.(all_λs)

        if ! bezier_radiative_transfer
            α5[i] = total_continuum_absorption([c_cgs/5e-5], layer.temp, nₑ, n_dict, partition_funcs)[1]
        end

        nₑ, n_dict, α_cntm_layer
    end
    nₑs = first.(triples)
    #put number densities in a dict of vectors, rather than a vector of dicts.
    n_dicts = getindex.(triples, 2)
    number_densities = Dict([spec=>[n[spec] for n in n_dicts] for spec in keys(n_dicts[1]) 
                             if spec != species"H III"])
    #vector of continuum-absorption interpolators
    α_cntm = last.(triples) 

    source_fn = blackbody.((l->l.temp).(atm.layers), all_λs')
    cntm = nothing
    if return_cntm
        cntm, _ = if bezier_radiative_transfer
            RadiativeTransfer.BezierTransfer.radiative_transfer(atm, α, source_fn, n_mu_points)
        else
            RadiativeTransfer.MoogStyleTransfer.radiative_transfer(atm, α, source_fn, α5, n_mu_points)
        end
    end

    if hydrogen_lines
        for (i, (layer, n_dict, nₑ)) in enumerate(zip(atm.layers, n_dicts, nₑs))
            hydrogen_line_absorption!(view(α, i, :), wl_ranges, layer.temp, nₑ,
                                      n_dict[species"H_I"],  n_dict[species"He I"],
                                      partition_funcs[species"H_I"](log(layer.temp)), vmic*1e5, 
                                      hydrogen_line_window_size*1e-8; 
                                      use_MHD=use_MHD_for_hydrogen_lines)
        end
    end

    line_absorption!(α, linelist, wl_ranges, [layer.temp for layer in atm.layers], nₑs,
        number_densities, partition_funcs, vmic*1e5, α_cntm, cutoff_threshold=line_cutoff_threshold)
    
    flux, intensity = if bezier_radiative_transfer
        RadiativeTransfer.BezierTransfer.radiative_transfer(atm, α, source_fn, n_mu_points)
    else
        RadiativeTransfer.MoogStyleTransfer.radiative_transfer(atm, α, source_fn, α5, n_mu_points)
    end

    # collect the indices corresponding to each wavelength range
    wl_lb_ind = 1 # the index into α of the lowest λ in the current wavelength range
    subspectra = []
    for λs in wl_ranges
        wl_inds = wl_lb_ind : wl_lb_ind + length(λs) - 1
        push!(subspectra, wl_inds)
        wl_lb_ind += length(λs)
    end

    # contruct a dict with strings instead of Korg.Species as keys (facilitates use from python)
    ns = Dict([string(spec) => n for (spec, n) in pairs(number_densities)])

    (flux=flux, cntm=cntm, intensity=intensity, alpha=α, number_densities=ns, 
    electron_number_density=nₑs, wavelengths=all_λs.*1e8, subspectra=subspectra)
end

"""
    format_A_X(default_metals_H, default_alpha_H, abundances; kwargs... )

Returns a 92 element vector containing abundances in ``A(X)`` (``\\log_{10}(X/H) + 12``) format for
elements from hydrogen to uranium.

# Arguments
You can specify abundance with these positional arguments.  All are optional, but if 
`default_alpha_H` is provided, `default_metals_H` must be as well. 
- `default_metals_H` (default: 0), i.e. [metals/H] is the ``\\log_{10}`` solar-relative abundance of elements heavier 
   than He. It is overriden by `default_alpha` and `abundances` on a per-element basis.  
- `default_alpha_H` (default: same as `default_metals_H`), i.e. [alpha/H] is the ``\\log_{10}`` 
   solar-relative abundance of the alpha elements (defined to be C, O, Ne, Mg, Si, S, Ar, Ca, and 
   Ti). It is overriden by `abundances` on a per-element basis.
- `abundances` is a `Dict` mapping atomic numbers or symbols to [``X``/H] abundances.  (Set 
  `solar_relative=false` to use ``A(X)`` abundances instead.) These override `default_metals_H`.
  This is the only way to specify an abundance of He that is non-solar.

# Keyword arguments
- `solar_relative` (default: true): When true, interpret abundances as being in \\[``X``/H\\] 
  (``\\log_{10}`` solar-relative) format.  When false, interpret them as ``A(X)`` abundances, i.e. 
   ``A(x) = \\log_{10}(n_X/n_\\mathrm{H}) + 12``, where ``n_X`` is the number density of ``X``.
   Note that abundances not specified default to the solar value still depend on the solar value, as
   they are set according to `default_metals_H` and `default_alpha_H`.
- `solar_abundances` (default: `Korg.asplund_2020_solar_abundances`) is the set of solar abundances to 
  use, as a vector indexed by atomic number. `Korg.asplund_2009_solar_abundances` and 
  `Korg.grevesse_2007_solar_abundances` are also provided for convienience.
"""
function format_A_X(default_metals_H::R1=0.0, default_alpha_H::R2=default_metals_H, 
                    abundances::Dict{K, V}=Dict{UInt8, Float64}();  
                    solar_relative=true, solar_abundances=default_solar_abundances
                    ) where {K, V, R1 <: Real, R2 <: Real}
    # make sure the keys of abundances are valid, and convert them to Z if they are strings
    clean_abundances = Dict{UInt8, V}()
    for (el, abund) in abundances
        if el isa AbstractString
            if ! (el in keys(Korg.atomic_numbers))
                throw(ArgumentError("$el isn't a valid atomic symbol."))
            elseif Korg.atomic_numbers[el] in keys(abundances)
                throw(ArgumentError("The abundances of $el was specified by both atomic number and atomic symbol."))
            else
                clean_abundances[Korg.atomic_numbers[el]] =  abund
            end
        elseif el isa Integer
            if ! (1 <= el <= MAX_ATOMIC_NUMBER)
                throw(ArgumentError("Z = $el is not a supported atomic number."))
            else
                clean_abundances[el] = abund
            end
        else
            throw(ArgumentError("$el isn't a valid element. Keys of the abundances dict should be strings or integers."))
        end
    end

    correct_H_abund = solar_relative ? 0.0 : 12.0
    if 1 in keys(clean_abundances) && clean_abundances[1] != correct_H_abund
        silly_abundance, silly_value = solar_relative ? ("[H/H]", 0) : ("A(H)", 12)
        throw(ArgumentError("$silly_abundance set, but $silly_abundance = $silly_value by " *
                            "definition. Adjust \"metallicity\" and \"abundances\" to implicitly " *
                            "set the amount of H"))
    end

    #populate A(X) vector
    alpha_els = [6, 8, 10, 12, 14, 16, 18, 20, 22]
    map(1:MAX_ATOMIC_NUMBER) do Z
        if Z == 1 #handle hydrogen
            12.0
        elseif Z in keys(clean_abundances) #if explicitely set
            if solar_relative
                clean_abundances[Z] + solar_abundances[Z]
            else
                clean_abundances[Z]
            end
        elseif Z in alpha_els
            solar_abundances[Z] + default_alpha_H
        else #if not set, use solar value adjusted for metallicity
            Δ = default_metals_H * (Z >= 3) #only adjust for metals, not H or He
            solar_abundances[Z] + Δ
        end
    end
end
# handle case where metallicity and alpha aren't specified but individual abundances are
format_A_X(abundances::Dict; kwargs...) = format_A_X(0, abundances; kwargs...)
# handle case where alpha isn't specified but individual abundances are
format_A_X(default_metallicity::R, abundances::Dict; kwargs...) where R <: Real = 
    format_A_X(default_metallicity, default_metallicity, abundances; kwargs...) 

"""
    get_metals_H(A_X)

Calculate [metals/H] given a vector, `A_X` of absolute abundances, ``A(X) = \\log_{10}(n_M/n_\\mathrm{H})``.
See also [`get_alpha_H`](@ref).
"""
function get_metals_H(A_X; solar_abundances=default_solar_abundances)
   _get_multi_X_H(A_X, 3:MAX_ATOMIC_NUMBER, solar_abundances)
end

"""
    get_alpha_H(A_X)

Calculate [α/H] given a vector, `A_X` of absolute abundances, ``A(X) = \\log_{10}(n_α/n_\\mathrm{H})``.
Here, the alpha elements are defined to be O, Ne, Mg, Si, S, Ar, Ca, Ti.  See also 
[`get_alpha_H`](@ref).
"""
function get_alpha_H(A_X; solar_abundances=default_solar_abundances)
    _get_multi_X_H(A_X, 8:2:22, solar_abundances)
end

"""
Given a vector of abundances, `A_X`, get [I+J+K/H], where `Zs = [I,J,K]` is a vector of atomic 
numbers.  This is used to calculate, for example, [α/H] and [metals/H].
"""
function _get_multi_X_H(A_X, Zs, solar_abundances)
    # there is no logsumexp in the julia stdlib, but it would make this more stable.
    # these are missing "+ 12", but it cancels out
    A_mX = log10(sum(10^A_X[Z] - 12 for Z in Zs))
    A_mX_solar = log10(sum(10^solar_abundances[Z] - 12 for Z in Zs))
    A_mX - A_mX_solar
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
