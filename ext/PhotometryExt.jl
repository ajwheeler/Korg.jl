module PhotometryExt

public compute_magnitude

import Interpolations: linear_interpolation
using Trapz
import PhotometricFilters

"""
    compute_magnitude(spec, photometric_filter)
    compute_magnitude(flux, wavelengths, photometric_filter)
    compute_magnitude(spec, filter_trans, filter_wave)
    compute_magnitude(flux, wavelengths, filter_trans, filter_wave)

Compute a photometric magnitude from a synthetic spectrum. This should only be used to create colors or magnitude differences (e.g. relative to Vega).

# Arguments
- `spec`, an object returned by `Korg.synthesize`. You can pass `flux` and `wavelengths` (in Å) directly instead (which are in Korg's default units).
- `photometric_filter`, a photometric filter object from `PhotometricFilters.jl`. You can pass `filter_trans` (the throughput) and `filter_wave` (in Å) directly instead.

Adapted from Eddie Schlafly.
"""
function compute_magnitude(spec, photometric_filter::PhotometricFilters.PhotometricFilter)
    compute_magnitude_core(spec.flux, spec.wavelengths, photometric_filter.throughput, photometric_filter.wave)
end
function compute_magnitude(spec_flux, spec_wave, photometric_filter::PhotometricFilters.PhotometricFilter)
    compute_magnitude_core(spec_flux, spec_wave, photometric_filter.throughput, photometric_filter.wave)
end
function compute_magnitude(spec, filter_trans, filter_wave)
    compute_magnitude_core(spec.flux, spec.wavelengths, filter_trans, filter_wave)
end
function compute_magnitude(spec_flux, spec_wave, filter_trans, filter_wave)
    compute_magnitude_core(spec_flux, spec_wave, filter_trans, filter_wave)
end

"""
Compute a magnitude in a given filter from a spectrum and filter transmission curve.
"""
function compute_magnitude_core(spec_flux, spec_wave, filter_trans, filter_wave)
    # add atmospheric transmission separately later
    # maybe I should pass a default resolution on which to compute the integral?
    @assert length(spec_wave) == length(spec_flux)

    # check spec_wave extends beyond filter_wave
    if minimum(filter_wave) < minimum(spec_wave) || maximum(filter_wave) > maximum(spec_wave)
        error("Spectrum does not cover the filter range")
    end
    # check spec_wave is longer than filter_wave
    if length(spec_wave) < length(filter_wave)
        throw(ArgumentError("Spectrum (length = $(length(spec_wave))) is lower resolution than the filter (length "*
            "= $(length(filter_wave)). Please use a higher resolution spectrum."))
    end

    # interpolate filter_wave to spec_wave
    interp_linear_extrap = linear_interpolation(filter_wave, filter_trans; extrapolation_bc=0)
    filter_trans_interp = interp_linear_extrap(spec_wave)

    # it is not clear to me that the units work out correctly here even if
    # one uses erg/s/cm^2/Å (which requires using get_radius)
    h = 6.62606957e-27  # erg * s
    c = 2.99792458e10  # cm/s
    flux_to_number = h * c ./ (spec_wave * 1e-8)  # erg

    spec_ref = 3631e-23 * c ./ (spec_wave * 1e-8) ./ spec_wave # working Mgy
    Iref = trapz(spec_wave, spec_ref .* filter_trans_interp .* flux_to_number)
    Ispec = trapz(spec_wave, spec_flux .* filter_trans_interp .* flux_to_number)

    return -2.5 * log10(Ispec / Iref) # AB magnitudes
end

# returns radius in cm given logg base 10 of cm/s^2
# assumes 1 solar mass
# conceptually use this as compute_magnitude(sol.flux.*(get_radius(4.43775056282).^2)./(1e8),sol.wavelengths,t.throughput,t.wave)
function get_radius(logg)
    # Constants
    G = 6.67430e-8  # gravitational constant in cgs
    M_sun = 1.989e33  # solar mass in g
    # Surface gravity in cgs from log(g)
    g = 10^logg # cm/s^2
    # Using g = GM/R², solve for R
    # R = sqrt(GM/g)
    radius = sqrt((G * M_sun) / g)  # in cm
    return radius
end

end # module
