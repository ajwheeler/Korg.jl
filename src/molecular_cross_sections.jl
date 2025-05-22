using Interpolations: interpolate!, Gridded, Linear, extrapolate
using HDF5

struct MolecularCrossSection
    wls::Wavelengths
    itp
    species::Korg.Species
end

"""
    MolecularCrossSection(linelist, wl_params...; cutoff_alpha=1e-30, log_temp_vals=3:0.025:5, verbose=true)

Precompute the molecular absorption cross section for a given linelist and set of wavelengths. The
`MolecularCrossSection` object can be passed to [`synthesize`](@ref) and potentially speed up the
calculation significantly.  At present, Korg only supports precomputed cross-sections created by
this function, though they can be saved and loaded using [`save_molecular_cross_section`](@ref) and
[`read_molecular_cross_section`](@ref).

# Arguments

  - `linelist`: A vector of `Line` objects representing the molecular linelist.  These must be of the
    same species.
  - `wl_params...`: Parameters specifying the wavelengths in the same format that [`synthesize`](@ref)
    expects.

# Keyword Arguments

  - `cutoff_alpha` (default: 1e-30): The value of the single-line absorption coefficient (in cm^-1) at
    which to truncate the profile.
  - `log_temp_vals` (default: 3:0.025:5): The log10 of the temperatures at which to precompute the
    cross-section.
  - `verbose` (default: true): Whether to print progress information.

!!! tip

    The default values of `vmic_vals`, `log_temp_vals`, and ` `cutoff_alpha` were chosen to ensure that lines in the APOGEE linelist ([`get_APOGEE_DR17_linelist`](@ref)) could be accurately
    reproduced (better than 10^-3 everywhere). You should verify that they yield acceptable accuracy
    for other applications by comparing spectra synthesize with and without precomputing the
    molecular cross-section.
"""
function MolecularCrossSection(linelist, wl_params...; cutoff_alpha=1e-32,
                               vmic_vals=[(0.0:1/3:1.0)...; 1.5; (2:2/3:(5+1/3))...],
                               log_temp_vals=3:0.04:5, verbose=true)
    wls = Wavelengths(wl_params...)
    all_specs = [l.species for l in linelist]
    if !all(Ref(all_specs[1]) .== all_specs)
        throw(ArgumentError("All lines must be of the same species"))
    end
    species = all_specs[1]

    α = zeros(length(vmic_vals), length(log_temp_vals), length(wls))

    # set both the continuum absorption coef (cntm) and the cutoff absorption coef to
    # unity.  Handle the cutoff value by scaling the number density of the molecule
    # (in n_dict).
    Ts = 10 .^ log_temp_vals
    nₑ = zeros(length(log_temp_vals))
    n_dict = Dict(species => 1 / cutoff_alpha)
    ξ = 0.0
    cntm = fill(λ -> 1.0, length(log_temp_vals))

    for (i, vmic) in enumerate(vmic_vals)
        ξ = vmic * 1e5 #km/s to cm/s
        Korg.line_absorption!(view(α, i, :, :), linelist, wls, Ts, nₑ, n_dict,
                              Korg.default_partition_funcs, ξ, cntm;
                              verbose=verbose, cutoff_threshold=1.0)
    end

    species = all_specs[1]
    itp = extrapolate(interpolate!((vmic_vals, log_temp_vals, wls), α .* cutoff_alpha,
                                   (Gridded(Linear()), Gridded(Linear()), Gridded(Linear()))), 0.0)
    MolecularCrossSection(wls, itp, species)
end

function Base.show(io::IO, sigma::MolecularCrossSection)
    Tmin, Tmax = round.(Int, 10 .^ sigma.itp.itp.knots[2][[1, end]])
    vmicmin, vmicmax = round.(Int, sigma.itp.itp.knots[1][[1, end]])
    wlmin, wlmax = round.(Int, sigma.wls[[1, end]] .* 1e8)
    println(io, "Molecular cross-section for ", sigma.species,
            " ($(wlmin) Å ≤ λ ≤ $(wlmax) Å, $(Tmin) K ≤ T ≤ $(Tmax) K, $(vmicmin) km/s ≤ vmic ≤ $(vmicmax) km/s)")
end

"""
    interpolate_molecular_cross_sections!(α, molecular_cross_sections, λs, Ts, vmic, number_densities)

Interpolate the molecular cross-sections and add them to the total absorption coefficient `α`.
See [`MolecularCrossSection`](@ref) for more information.
"""
function interpolate_molecular_cross_sections!(α::AbstractArray{R}, molecular_cross_sections,
                                               λs::Wavelengths, Ts, vmic,
                                               number_densities) where R<:Real
    if length(molecular_cross_sections) == 0
        return
    end

    for sigma in molecular_cross_sections
        # this is an inefficient order in which to write to α, but doing the interpolations this
        # way uses less memory.
        for i in 1:size(α, 1)
            vm = vmic isa Number ? vmic : vmic[i]
            α[i, :] .+= sigma.itp.(vm, log10(Ts[i]), λs) * number_densities[sigma.species][i]
        end
    end
end

"""
    save_molecular_cross_section(filename, cross_section)

Save a precomputed molecular cross-section to a file.
See also [`MolecularCrossSection`](@ref), [`read_molecular_cross_section`](@ref).
"""
function save_molecular_cross_section(filename, cross_section)
    wls = cross_section.wls
    itp = cross_section.itp
    species = cross_section.species

    HDF5.h5open(filename, "w") do file
        HDF5.write(file, "wls", [(l[begin], step(l), l[end]) for l in wls.wl_ranges])
        HDF5.write(file, "vmic_vals", collect(itp.itp.knots[1]))
        HDF5.write(file, "T_vals", collect(itp.itp.knots[2]))
        HDF5.write(file, "vals", itp.itp.coefs)
        HDF5.write(file, "species", string(species))
    end
end

"""
    read_molecular_cross_section(filename)

Read a precomputed molecular cross-section from a file created by
[`save_molecular_cross_section`](@ref).
See also [`MolecularCrossSection`](@ref).
"""
function read_molecular_cross_section(filename)
    HDF5.h5open(filename, "r") do file
        wls = map(HDF5.read(file, "wls")) do (start, step, stop)
            start:step:stop
        end |> Wavelengths
        vmic_vals = HDF5.read(file, "vmic_vals")
        logT_vals = HDF5.read(file, "T_vals")
        # doesn't actually matter that it's mmaped because it's passed to interpolate!
        alpha_vals = HDF5.readmmap(file["vals"])
        species = Species(HDF5.read(file, "species"))

        itp = extrapolate(interpolate!((vmic_vals, logT_vals, wls), alpha_vals,
                                       (Gridded(Linear()), Gridded(Linear()), Gridded(Linear()))),
                          0.0)

        MolecularCrossSection(wls, itp, species)
    end
end
