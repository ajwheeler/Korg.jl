
import Interpolations: linear_interpolation 
using Trapz, DelimitedFiles

# add atmospheric transmission separately later
# maybe I should pass a default resolution on which to compute the integral?
# adapted from Eddie Schlafly
function compute_magnitude(spec_flux, spec_wave, filter_trans, filter_wave)
    # check spec_wave extends beyond filter_wave
    if minimum(filter_wave) < minimum(spec_wave) || maximum(filter_wave) > maximum(spec_wave)
        error("Spectrum does not cover the filter range")
    end
    # check spec_wave is longer than filter_wave
    if length(spec_wave) < length(filter_wave)
        error("Spectrum is lower resolution than the filter")
    end

    # interpolate filter_wave to spec_wave
    interp_linear_extrap = linear_interpolation(filter_wave, filter_trans, extrapolation_bc=0) 
    filter_trans_interp = interp_linear_extrap(spec_wave)

    # it is not clear to me that the units work out correctly here even if
    # one uses erg/s/cm^2/Å (which requires using get_radius)
    h = 6.62606957e-27  # erg * s
    c = 2.99792458e10  # cm/s
    flux_to_number = h*c ./(spec_wave*1e-8);  # erg

    spec_ref = 3631e-23*c ./(spec_wave*1e-8)./spec_wave # working Mgy
    Iref = trapz(spec_wave,spec_ref.*filter_trans_interp.*flux_to_number)
    Ispec = trapz(spec_wave,spec_flux.*filter_trans_interp.*flux_to_number)
    
    return -2.5 * log10(Ispec/Iref) # AB magnitudes
end

# returns radius in cm given logg base 10 of cm/s^2
# assumes 1 solar mass
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

# takes filter name and returns lambda and throughput (including atm transparancy)
function parse_DECam_filter(filter_name)
    base = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/2024_11_18/Korg.jl/data/filter_curves/DECam/DECam_"
    data = readdlm(base * filter_name * ".txt"; comments=true)
    # msk rows where data[:,2] is zero (no transmission)
    msk = data[:, 2] .!= 0
    return data[msk, 1], data[msk, 2]
end