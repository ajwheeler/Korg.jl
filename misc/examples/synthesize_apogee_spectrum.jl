using Korg

# this file was created by creating a Korg.MolecularCrossSection
# and saved with Korg.save_molecular_cross_section
water = Korg.read_molecular_cross_section("apogee_water_sigma.h5")

# stellar params
Teff = 5777
logg = 4.44
vmic = 1.2

# abundances
M_H = -0.1
alpha_H = 0.2
C_H = -0.3

# constract length-92 vector of A(X)-format abundances
A_X = format_A_X(M_H, alpha_H, Dict("C" => C_H))

# create model atmosphere by interpolating from one of a few grid
# (RGB will be normal SDSS MARCS grid)
atm = interpolate_marcs(Teff, logg, A_X)

# set up sparse matrix that applies the LSF and resamples to the apogee grid
synthesis_wavelengths = 15_000:0.01:17_000
apogee_wavelengths = 10 .^ range(; start=log10(15100.802), step=6e-6, length=8575)
LSF = Korg.compute_LSF_matrix(synthesis_wavelengths, apogee_wavelengths, 22_500)

# don't include water lines because we will use the precomputed cross-section instead
linelist = Korg.get_APOGEE_DR17_linelist(; include_water=false)

# this returns a "solution" object with lots of data
# @time is a macro that reports how long this line takes
@time sol = synthesize(atm, linelist, A_X, synthesis_wavelengths;
                       vmic=vmic, use_MHD_for_hydrogen_lines=false,
                       molecular_cross_sections=[water])

# get continuum normalized flux by taking ratio of flux and theoretical continuum
# (synth does this for you)
# use sol.flux if you want to work with unnormalized spectra
rectified_flux = sol.flux ./ sol.cntm

final_spectrum = LSF * rectified_flux
