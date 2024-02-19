using Interpolations: GriddedInterpolation, interpolate, Gridded, NoInterp, Linear
using DSP: gaussian

struct MolecularCrossSection
    itp :: GriddedInterpolation
    species :: Korg.Species

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
    function MolecularCrossSection(linelist, wls; 
                                   cutoff_alpha=1e-30, log_temp_vals=3:0.025:5, verbose=true)
        all_specs = [l.species for l in linelist]
        if !all(Ref(all_specs[1]) .== all_specs)
            throw(ArgumentError("All lines must be of the same species"))
        end

        α = zeros(length(log_temp_vals), sum(length.(wls)))
        
        species = all_specs[1]

        # set both the continuum absorption coef (cntm) and the cutoff absorption coef to 
        # unity.  Handle the cutoff value by scaling the number density of the molecule
        # (in n_dict).
        Ts = 10 .^ log_temp_vals
        nₑ = zeros(length(log_temp_vals))
        n_dict = Dict(species => 1/cutoff_alpha)
        ξ = 0.0
        cntm = fill(λ -> 1.0, length(log_temp_vals))
        
        Korg.line_absorption!(α, linelist, wls*1e-8, Ts, nₑ, n_dict, Korg.default_partition_funcs, ξ, cntm; 
                              verbose=verbose, cutoff_threshold=1.0)

        species = all_specs[1]
        itp = interpolate((log_temp_vals, 1:size(α, 2)),
                          α .* cutoff_alpha, 
                          (Gridded(Linear()), NoInterp()))
        new(itp, species)
    end
end

"""
"""
function interpolate_molecular_cross_sections!(α::Matrix{R}, molecular_cross_sections, Ts, vmic,
                                               wls, number_densities) where R <: Real
    if length(molecular_cross_sections) == 0
        return
    end

    α_mol = zeros(R, size(α))
    for sigma in molecular_cross_sections
        α_mol .+= sigma.itp.(log10.(Ts), (1:size(α, 2))') .* number_densities[sigma.species]
    end
    α .+= α_mol

    if vmic != 0
        throw(ArgumentError("nonzero vmic not yet supported for molecular cross sections"))
    #    for wl_range in wls

    #    end
    end
    ;
end

#function apply_smoothing_filter!(flux, kernel_width=150;
#    kernel_pixels=301)
#    
#    @assert isodd(kernel_pixels) # required for offset math to work
#    sampled_kernel = gaussian(kernel_pixels, kernel_width/kernel_pixels)
#    sampled_kernel ./= sum(sampled_kernel) # normalize (everything should cancel in the end, but I'm paranoid)
#    offset = (length(sampled_kernel) - 1) ÷ 2
#
#    map(region_inds) do r
#        buffered_region = r[1]-offset : r[end]+offset
#        flux[buffered_region] = conv(sampled_kernel, flux[buffered_region])[offset+1 : end-offset]
#    end
#end