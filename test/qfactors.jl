@testset "qfactors" begin
    # Check a Q factor gives the right RV precision
    @test abs(Korg.precRV_fromQ(1100,100,8000)-30.470741601352298) < 1e-10

    # Check that a solar spectrum at APOGEE resolution gives the expected Q factor
    delLog = 6e-6; 
    wavetarg = 10 .^range((start=4.179-125*delLog),step=delLog,length=8575+125)
    minw, maxw = extrema(wavetarg);
    
    wl_lo, wl_hi, wl_step = 15_000, 17_000, 1//100
    x_model = wl_lo:wl_step:wl_hi

    LSF_model = Korg.compute_LSF_matrix(x_model, wavetarg, 22_500, renormalize_edge=true, verbose=false)

    apolines = Korg.get_APOGEE_DR17_linelist(include_water=true);

    A_X = copy(Korg.grevesse_2007_solar_abundances)
    atm = Korg.interpolate_marcs(5777, 4.4, Korg.grevesse_2007_solar_abundances)
    kout = synthesize(atm, apolines, A_X, wl_lo, wl_hi, wl_step, vmic=0, hydrogen_lines=true, hydrogen_line_window_size=300, electron_number_density_warn_threshold=1e10);
    tspec = (kout.flux)./(kout.cntm);
    @test abs(Korg.Qfactor(tspec, x_model, wavetarg, LSF_model)-1182.420317367743)<1e-10

    msk = ones(Bool,length(wavetarg))
    msk[1:100].=false
    msk[end-100:100].=true
    RV_prec = Korg.precRV_MeasuredNoise(tspec, x_model, wavetarg, LSF_model, 1 ./100*ones(length(wavetarg)),msk_lres=msk)
    @test abs(RV_prec-29.490659698846365)<1e-10
    @test abs(Korg.precRV_fromQ(Q,100,count(msk)).-27.34006821159959)<1e-10
    # I think these two should agree... but the problem must be from the continuum defintion and how that propagates to S/N
end