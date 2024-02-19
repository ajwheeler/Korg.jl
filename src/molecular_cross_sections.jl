using Interpolations: GriddedInterpolation, interpolate, Gridded, NoInterp, Linear
using DSP: gaussian
using HDF5

struct MolecularCrossSection
    wls # inexpensive to save, just for debugging
    itp :: GriddedInterpolation
    species :: Korg.Species
end

"""
    MolecularCrossSection(linelist, wls; cutoff_alpha=1e-30, log_temp_vals=3:0.025:5, verbose=true)

Precompute the molecular absorption cross section for a given linelist and set of wavelengths. The 
resulting object can be passed to [`synthesize`](@ref) and potentially speed up the calculation 
significantly.  At present, Korg only supports opacity tables computed by this function.

# Arguments
- `linelist`: A vector of `Line` objects representing the molecular linelist.  These 
    must be of the same species.
- `wls`: A vector of wavelength ranges (in Å) at which to precompute the cross section.  *These must
  match the wavelengths used for any subsequent synthesis exactly*.

# Keyword Arguments
- `cutoff_alpha` (default: 1e-30): The value of the single-line absorption coefficient (in cm^-1) at 
   which to truncate the profile.
- `log_temp_vals` (default: 3:0.025:5): The log10 of the temperatures at which to precompute the 
   cross-section.
- `verbose` (default: true): Whether to print progress information.

!!! tip
    The default values of `log_temp_vals` and `cutoff_alpha` were chosen by ensuring that water 
    lines in the APOGEE linelist ([`get_APOGEE_DR17_linelist`](@ref)) could be accurately reproduced 
    (better than 10^-3 everywhere). You should verify that they yield acceptable accuracy for other 
    applications by comparing spectra synthesize with and without precomputing the molecular 
    cross-section.

# Returns
An interpolator object which can be passed to [`synthesize`](@ref).
"""
function MolecularCrossSection(linelist, wls; cutoff_alpha=1e-32, 
                               vmic_vals=[(0.0:1/3:1.0)... ; 1.5 ; (2:2/3:(5+1/3))...],
                               log_temp_vals=3:0.04:5, verbose=true)
    all_specs = [l.species for l in linelist]
    if !all(Ref(all_specs[1]) .== all_specs)
        throw(ArgumentError("All lines must be of the same species"))
    end

    α = zeros(length(vmic_vals), length(log_temp_vals), sum(length.(wls)))
        
    species = all_specs[1]

    # set both the continuum absorption coef (cntm) and the cutoff absorption coef to 
    # unity.  Handle the cutoff value by scaling the number density of the molecule
    # (in n_dict).
    Ts = 10 .^ log_temp_vals
    nₑ = zeros(length(log_temp_vals))
    n_dict = Dict(species => 1/cutoff_alpha)
    ξ = 0.0
    cntm = fill(λ -> 1.0, length(log_temp_vals))
        
    for (i, vmic) in enumerate(vmic_vals)
        ξ = vmic * 1e5 #km/s to cm/s
        Korg.line_absorption!(view(α, i, :, :), linelist, wls*1e-8, Ts, nₑ, n_dict,
                              Korg.default_partition_funcs, ξ, cntm; 
                              verbose=verbose, cutoff_threshold=1.0)
    end

    species = all_specs[1]
    itp = interpolate((vmic_vals, log_temp_vals, 1:size(α, 3)),
                      α .* cutoff_alpha, 
                      (Gridded(Linear()), Gridded(Linear()), NoInterp()))
    MolecularCrossSection(wls, itp, species)
end


"""
    interpolate_molecular_cross_sections!(α, molecular_cross_sections, Ts, vmic, number_densities)

Interpolate the molecular cross-sections and add them to the total absorption coefficient `α`.
See [`MolecularCrossSection`](@ref) for more information.
"""
function interpolate_molecular_cross_sections!(α::Matrix{R}, molecular_cross_sections, Ts, vmic,
                                               number_densities) where R <: Real
    if length(molecular_cross_sections) == 0
        return
    end

    α_mol = zeros(R, size(α))
    for sigma in molecular_cross_sections
        α_mol .+= sigma.itp.(vmic, log10.(Ts), (1:size(α, 2))') .* number_densities[sigma.species]
    end
    α .+= α_mol
    ;
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
        HDF5.write(file, "wls", [(l[begin], step(l), l[end]) for l in wls])
        HDF5.write(file, "vmic_vals", collect(itp.knots[1]))
        HDF5.write(file, "T_vals", collect(itp.knots[2]))
        HDF5.write(file, "vals", itp.coefs)
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
        end
        vmic_vals = HDF5.read(file, "vmic_vals")
        T_vals = HDF5.read(file, "T_vals")
        alpha_vals = HDF5.read(file, "vals")
        species = Species(HDF5.read(file, "species"))

        itp = interpolate((vmic_vals, T_vals, 1:sum(length.(wls))),
                          alpha_vals,
                          (Gridded(Linear()), Gridded(Linear()), NoInterp()))

        MolecularCrossSection(wls, itp, species)
    end
end