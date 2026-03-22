# This notebook demonstrates how to compute theoretical bounds on RV-precision given a synthetic
# spectrum and an LSF.  Korg has a few functions that make this easier, based on [Bouchy et al.
# 2001, A&A, 374, 733](https://ui.adsabs.harvard.edu/abs/2001A%26A...374..733B/abstract).

# In addition to Korg, we'll use `PythonPlot` in this example. If you don't have it installed, you
# can run `using Pkg; Pkg.add("PythonPlot")` to install it.

using Suppressor # hide #md
@suppress begin # hide #md
    using Korg, PythonPlot
end # hide #md

# First, specify the wavelengths and line spread function (LSF).
# We'll use the APOGEE instument's wavelength coverage and resolution for this example.
# Reduced APOGEE spectra are resampled onto wavelengths which are uniform in log-wavelength.

delLog = 6e-6;
apowls = 10 .^ range((start = 4.179 - 125 * delLog); step=delLog, length=8575 + 125)

# We compute the LSF from the APOGEE resolution (``R \approx 22,500``) as a sparse matrix.
# (See the [LSF explanation here](@ref LSF) for more details.)

LSF = Korg.compute_LSF_matrix((15_000, 17_000), apowls, 22_500)

# Now, let's synthesize a spectrum, and simultaneously apply the LSF and resample to the APOGEE
# wavelength grid.

apolines = Korg.get_APOGEE_DR17_linelist(; include_water=true)
wls, flux, _ = Korg.synth(; linelist=apolines,
                          wavelengths=(15_000, 17_000),
                          Teff=5777,
                          logg=4.44,
                          M_H=-1.1,)

SNR = 50.0
obs_err = (LSF * flux) ./ SNR;

# # Exact radial velocity precision
#
# This function provides the exact RV precision bound (in the small redshift limit) in m/s. Note
# that we are using the unconvolved flux here (`flux`), not the values with the LSF applied
# (`LSF * flux`).

noise_prec = Korg.RV_prec_from_noise(flux, (15_000, 17_000), apowls, LSF, obs_err)

# # ``Q``-factor calculation
# We can also calculate the ``Q``-factor, which provides an approximately ``S/N``-independent notion
# of the radial velocity "information" in the spectrum.

Q = Korg.Qfactor(flux, (15_000, 17_000), apowls, LSF)

# We can calculate the appriximate RV precision bound from ``Q``.  Technically,
# the precision depends on the root-mean-squared per-pixel ``S/N`` value, but for non-line-blanket
# spectra, this is approximately equal to the ``S/N``. (Verifying this fact is left as an excercise
# to the reader.)

SNR_RMS = SNR # approximately
Npixels = length(apowls)
Q_prec = Korg.RV_prec_from_Q(Q, SNR, Npixels)

# This is a few m/s off of the exact calculation, but it's pretty close!

# # ``Q``-factor as a function of metallicity
#
# Let's examine how the ``Q``-factor changes as we change the metallicity for a star with solar
# ``T_\mathrm{eff}`` and ``\log g``. As we add more and deeper lines, there is more for an RV
# measurement to grab onto.

M_Hs = -5:1.0:0 #metallicity ([M/H]) values
Qs = map(M_Hs) do M_H
    _, flux, _ = Korg.synth(; linelist=apolines, wavelengths=(15_000, 17_000), Teff=5777, logg=4.44,
                            M_H=M_H,)
    Korg.Qfactor(flux, (15_000, 17_000), apowls, LSF)
end

figure(; figsize=(12, 4)) # hide #md
scatter(M_Hs, Qs)
xlabel("[M/H]")
ylabel(L"$Q$ [m/s]")
gcf() # hide #md

# We can use the ``Q``-factors to calculate the RV precision across ``S/N``. Note that these
# precision estimates are not the whole story at high ``S/N`` because effects like granulation become
# important.

figure(; figsize=(12, 4)) # hide #md
SNRs = 1:1:1000
for (M_H, Q) in zip(M_Hs, Qs)
    plot(SNRs, Korg.RV_prec_from_Q.(Q, SNRs, Npixels); label="[M/H] = $M_H")
end
legend()
ylabel("RV precision [m/s]")
xlabel("S/N")
yscale("log")
xscale("log")
gcf() # hide #md
