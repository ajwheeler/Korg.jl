using Interpolations: interpolate!, Gridded, Linear, extrapolate
using SparseArrays: rowvals, nzrange
using HDF5

struct MolecularCrossSection
    wls::Wavelengths
    itp
    species::Korg.Species
end

const _VMIC_GRID_DEPRECATION_MSG = "vmic-gridded molecular cross-sections are deprecated. " *
                                   "Microturbulence is now applied by convolving a vmic=0 cross-section at synthesis time, " *
                                   "which is faster to precompute, smaller to store, and more accurate for vmic values off " *
                                   "the precomputation grid.  Regenerate your cross-section with " *
                                   "MolecularCrossSection(linelist, wls) (without vmic_vals) to use the new format."

"""
    is_vmic_gridded(cross_section)

Whether a [`MolecularCrossSection`](@ref) is in the deprecated format which grids over vmic, rather
than the current format, which is computed at vmic=0 and broadened at synthesis time.
"""
is_vmic_gridded(σ::MolecularCrossSection) = length(σ.itp.itp.knots) == 3

"""
    MolecularCrossSection(linelist, wl_params...; kwargs...)

Precompute the molecular absorption cross section for a given linelist and set of wavelengths. The
`MolecularCrossSection` object can be passed to [`synthesize`](@ref) and potentially speed up the
calculation significantly.  At present, Korg only supports precomputed cross-sections created by
this function, though they can be saved and loaded using [`save_molecular_cross_section`](@ref) and
[`read_molecular_cross_section`](@ref).

# Arguments

  - `linelist`: A vector of `Line` objects representing the molecular linelist.  These must be of the
    same species.
  - `wl_param`: Parameters specifying the wavelengths in the same format that [`synthesize`](@ref)
    expects.  Because microturbulence is applied by convolution at synthesis time, the wavelength range of
    the cross-section should extend a few Å beyond the range you plan to synthesize.
    Otherwise absorption from just outside the table will be missing near its edges.

# Keyword Arguments

  - `cutoff_alpha` (default: 1e-32): The value of the single-line absorption coefficient (in cm^-1) at
    which to truncate the profile.
  - `log_temp_vals` (default: 3:0.04:5): The log10 of the temperatures at which to precompute the
    cross-section.
  - **deprecated** `vmic_vals` (default: `nothing`): If provided, a vector of microturbulence values
    (in km/s) at which to precompute the cross-section, producing a table in the deprecated
    vmic-gridded format (which is interpolated in vmic at synthesis time, rather than convolved).

!!! warning

    While this is part of the public API, it is not considered stable and may change in the future
    without a major version bump.

!!! tip

    The default values of `log_temp_vals` and `cutoff_alpha` were chosen to ensure that lines in
    the APOGEE linelist ([`get_APOGEE_DR17_linelist`](@ref)) could be accurately reproduced (better
    than 10^-3 everywhere). You should verify that they yield acceptable accuracy for other
    applications by comparing spectra synthesize with and without precomputing the molecular
    cross-section.
"""
function MolecularCrossSection(linelist, wl_params; cutoff_alpha=1e-32, vmic_vals=nothing,
                               log_temp_vals=3:0.04:5)
    wls = Wavelengths(wl_params)
    all_specs = [l.species for l in linelist]
    if !all(Ref(all_specs[1]) .== all_specs)
        throw(ArgumentError("All lines must be of the same species"))
    end
    species = all_specs[1]

    # set both the continuum absorption coef (cntm) and the cutoff absorption coef to
    # unity.  Handle the cutoff value by scaling the number density of the molecule
    # (in n_dict).
    Ts = 10 .^ log_temp_vals
    nₑ = zeros(length(log_temp_vals))
    # Korg will compute the effective broadener number density for vdW broadening, which requires
    # number density values for H I, He I, and H2. It will not be used, since only molecular lines
    # are present.
    n_dict = Dict(species => 1 / cutoff_alpha,
                  species"H I" => NaN, species"H2" => NaN, species"He I" => NaN)
    cntm = fill(λ -> 1.0, length(log_temp_vals))

    itp = if isnothing(vmic_vals)
        α = zeros(length(log_temp_vals), length(wls))
        Korg.line_absorption!(α, linelist, wls, Ts, nₑ, n_dict,
                              Korg.default_partition_funcs, 0.0, cntm; cutoff_threshold=1.0)
        extrapolate(interpolate!((log_temp_vals, wls), α .* cutoff_alpha,
                                 (Gridded(Linear()), Gridded(Linear()))), 0.0)
    else
        @warn _VMIC_GRID_DEPRECATION_MSG maxlog = 1
        α = zeros(length(vmic_vals), length(log_temp_vals), length(wls))
        for (i, vmic) in enumerate(vmic_vals)
            ξ = vmic * 1e5 #km/s to cm/s
            Korg.line_absorption!(view(α, i, :, :), linelist, wls, Ts, nₑ, n_dict,
                                  Korg.default_partition_funcs, ξ, cntm; cutoff_threshold=1.0)
        end
        extrapolate(interpolate!((vmic_vals, log_temp_vals, wls), α .* cutoff_alpha,
                                 (Gridded(Linear()), Gridded(Linear()), Gridded(Linear()))), 0.0)
    end
    MolecularCrossSection(wls, itp, species)
end

function Base.show(io::IO, sigma::MolecularCrossSection)
    wlmin, wlmax = round.(Int, sigma.wls[[1, end]] .* 1e8)
    if is_vmic_gridded(sigma)
        Tmin, Tmax = round.(Int, 10 .^ sigma.itp.itp.knots[2][[1, end]])
        vmicmin, vmicmax = round.(Int, sigma.itp.itp.knots[1][[1, end]])
        println(io, "Molecular cross-section for ", sigma.species,
                " ($(wlmin) Å ≤ λ ≤ $(wlmax) Å, $(Tmin) K ≤ T ≤ $(Tmax) K, " *
                "$(vmicmin) km/s ≤ vmic ≤ $(vmicmax) km/s, deprecated vmic-gridded format)")
    else
        Tmin, Tmax = round.(Int, 10 .^ sigma.itp.itp.knots[1][[1, end]])
        println(io, "Molecular cross-section for ", sigma.species,
                " ($(wlmin) Å ≤ λ ≤ $(wlmax) Å, $(Tmin) K ≤ T ≤ $(Tmax) K)")
    end
end

# the equivalent resolving power (R = λ/FWHM) of the microturbulent broadening kernel.
# Korg's doppler_width convention is σ_λ = λ₀ sqrt(kT/m + ξ²/2)/c, so broadening a vmic=0
# cross-section to microturbulence ξ requires a Gaussian with velocity width σ_v = ξ/√2.
_vmic_equivalent_R(ξ) = c_cgs / (2sqrt(log(2)) * ξ)

"""
    interpolate_molecular_cross_sections!(α, molecular_cross_sections, λs, Ts, vmic, number_densities;
                                          vmic_window_size=4.0)

Evaluate the molecular cross-sections at each wavelength and layer and add them to the total
absorption coefficient, `α`.

Microturbulent broadening is applied via (sparse) matrix multiplication (convolution in velocity
space), extending `vmic_window_size` standard deviations from the center.

See [`MolecularCrossSection`](@ref) for more information.
"""
function interpolate_molecular_cross_sections!(α::AbstractArray{R}, molecular_cross_sections,
                                               λs::Wavelengths, Ts, vmic,
                                               number_densities; vmic_window_size=4.0) where R<:Real
    if length(molecular_cross_sections) == 0
        return
    end
    n_layers = size(α, 1)

    for sigma in molecular_cross_sections
        n_species = number_densities[sigma.species]

        # deprecated, backwards-compatible branch
        if is_vmic_gridded(sigma)
            @warn _VMIC_GRID_DEPRECATION_MSG maxlog = 1
            # this is an inefficient order in which to write to α, but doing the interpolations
            # this way uses less memory.
            for i in 1:n_layers
                vm = vmic isa Number ? vmic : vmic[i]
                α[i, :] .+= sigma.itp.(vm, log10(Ts[i]), λs) * n_species[i]
            end
            continue
        end

        # one broadening matrix per unique vmic value, shared across layers
        unique_vmics = vmic isa Number ? [vmic] : unique(vmic)
        for vm in unique_vmics
            # at vmic == 0 the kernel is a delta function, so interpolate directly.
            # when vmic == 0, ∂α/∂vmic == 0 as well, so it's legit to drop any partials of vmic
            if vm == 0
                for i in 1:n_layers
                    (vmic isa Number ? vmic : vmic[i]) == vm || continue
                    α[i, :] .+= sigma.itp.(log10(Ts[i]), λs) * n_species[i]
                end
                continue
            end

            # Use the LSF machinery to compure a sparse matrix that applies a velocity-space 
            # gaussian convolution and resamples from the cross-section λs to the synthesis λs.
            vmic_matrix = compute_LSF_matrix(sigma.wls, λs, _vmic_equivalent_R(vm * 1e5);
                                   window_size=vmic_window_size, verbose=false)

            # keep track of which wavelengths get nothing from the matrix, becaue they are outside 
            # the precomputed cross-section's wavelength regions
            row_has_kernel = falses(length(λs))
            row_has_kernel[rowvals(vmic_matrix)] .= true
            fallback_inds = findall(!, row_has_kernel)

            # restrict the matrix product to the table columns under some kernel.  This makes
            # requests for a few wavelengths (α at the reference wavelength for anchored radiative
            # transfer) cheap: only a small piece of the table is interpolated per layer.
            cols = let nonempty = j -> !isempty(nzrange(vmic_matrix, j))
                first_col = findfirst(nonempty, 1:size(vmic_matrix, 2))
                isnothing(first_col) ? (1:0) : (first_col:findlast(nonempty, 1:size(vmic_matrix, 2)))
            end
            Msub = vmic_matrix[:, cols]
            sub_wls = sigma.wls[cols]

            for i in 1:n_layers
                (vmic isa Number ? vmic : vmic[i]) == vm || continue
                logT = log10(Ts[i])
                if !isempty(cols)
                    α[i, :] .+= (Msub * sigma.itp.(logT, sub_wls)) .* n_species[i]
                end
                for j in fallback_inds
                    α[i, j] += sigma.itp(logT, λs[j]) * n_species[i]
                end
            end
        end
    end
end

"""
    save_molecular_cross_section(filename, cross_section)

Save a precomputed molecular cross-section to an HDF5 file.
See also [`MolecularCrossSection`](@ref), [`read_molecular_cross_section`](@ref).

!!! warning

    While this is part of the public API, it is not considered stable and may change in the future
    without a major version bump.
"""
function save_molecular_cross_section(filename, cross_section)
    wls = cross_section.wls
    itp = cross_section.itp
    species = cross_section.species

    HDF5.h5open(filename, "w") do file
        HDF5.write(file, "wls", [(l[begin], step(l), l[end]) for l in wls.wl_ranges])
        if is_vmic_gridded(cross_section)
            @warn _VMIC_GRID_DEPRECATION_MSG maxlog = 1
            HDF5.write(file, "vmic_vals", collect(itp.itp.knots[1]))
            HDF5.write(file, "T_vals", collect(itp.itp.knots[2]))
        else
            HDF5.write(file, "T_vals", collect(itp.itp.knots[1]))
        end
        HDF5.write(file, "vals", itp.itp.coefs)
        HDF5.write(file, "species", string(species))
    end
end

"""
    read_molecular_cross_section(filename)

Read a precomputed molecular cross-section from a file created by
[`save_molecular_cross_section`](@ref).
See also [`MolecularCrossSection`](@ref).

!!! warning

    While this is part of the public API, it is not considered stable and may change in the future
    without a major version bump.
"""
function read_molecular_cross_section(filename)
    HDF5.h5open(filename, "r") do file
        wls = map(HDF5.read(file, "wls")) do (start, step, stop)
            # reconstruct with an explicit length so that floating-point noise in the stored
            # (cm-scale) endpoints can't add or drop a point
            n = round(Int, (stop - start) / step) + 1
            range(start * 1e8, stop * 1e8; length=n)
        end |> Wavelengths
        logT_vals = HDF5.read(file, "T_vals")
        # doesn't actually matter that it's mmaped because it's passed to interpolate!
        alpha_vals = HDF5.readmmap(file["vals"])
        species = Species(HDF5.read(file, "species"))

        itp = if haskey(file, "vmic_vals")
            @warn _VMIC_GRID_DEPRECATION_MSG maxlog = 1
            vmic_vals = HDF5.read(file, "vmic_vals")
            extrapolate(interpolate!((vmic_vals, logT_vals, wls), alpha_vals,
                                     (Gridded(Linear()), Gridded(Linear()), Gridded(Linear()))),
                        0.0)
        else
            extrapolate(interpolate!((logT_vals, wls), alpha_vals,
                                     (Gridded(Linear()), Gridded(Linear()))), 0.0)
        end

        MolecularCrossSection(wls, itp, species)
    end
end
