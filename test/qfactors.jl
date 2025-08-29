
@testset "qfactors" begin
    using Statistics: mean

    # Check a Q factor gives the right RV precision
    @test Korg.RV_prec_from_Q(1100, 100, 8000)≈30.470741601352298 atol=1e-10

    # wavelengths and LSF
    wl_lo, wl_hi = 15_100, 15_200

    delLog = 6e-6
    apowls = 10 .^ range((start = 4.179 - 125 * delLog); step=delLog, length=8575 + 125)
    apowls = apowls[wl_lo.<apowls.<wl_hi]

    LSF_model = Korg.compute_LSF_matrix((wl_lo, wl_hi), apowls, 22_500; verbose=false)

    # it's less accurate to not include water lines, but that's fine for this test
    apolines = Korg.get_APOGEE_DR17_linelist(; include_water=false)
    A_X = format_A_X()
    atm = Korg.interpolate_marcs(5777, 4.4, A_X)
    sol = synthesize(atm, apolines, A_X, (wl_lo, wl_hi))

    synth_flux = (sol.flux) ./ (sol.cntm)
    # this denominator is a normalization factor
    flux = (LSF_model * synth_flux)

    flux = (LSF_model * synth_flux) ./ sum(LSF_model; dims=2)[:]
    # for the Q factor calculation to be exact, the noise must be photon-dominated
    obs_err = sqrt.(flux) .* 0.01

    msk = ones(Bool, length(apowls))
    msk[1:100] .= false

    Q = Korg.Qfactor(synth_flux, (wl_lo, wl_hi), apowls, LSF_model; obs_mask=msk)
    @test Q≈877.6 atol=1
    SNR = flux ./ obs_err
    RMS_SNR = sqrt(mean(SNR[msk] .^ 2))
    Q_prec = Korg.RV_prec_from_Q(Q, RMS_SNR, count(msk))

    noise_prec = Korg.RV_prec_from_noise(synth_flux, (wl_lo, wl_hi), apowls, LSF_model, obs_err;
                                         obs_mask=msk)

    @test Q_prec ≈ noise_prec
end
