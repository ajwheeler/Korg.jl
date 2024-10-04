@testset "Wavelengths" begin
    wls = Korg.Wavelengths(15000, 15500)
    @test wls == Korg.Wavelengths(15000:0.01:15500)
    @test wls == Korg.Wavelengths([15000:0.01:15500])
    @test wls == Korg.Wavelengths(collect(15000:0.01:15500))

    # test automatic air to vacuum conversion
    @test Korg.Wavelengths(15000:0.01:15500; air_wavelengths=true).wl_range ==
          Korg.air_to_vacuum.(15000:0.01:15500)

    #TODO do these make sense to do?
    #@testset "synthesize integration test" begin
    #    # these are essentially tests of Wavelegnths, but we test "through"
    #    # synthesize because it's crucial that synthesize works

    #    # truncated atmosphere for testing
    #    atm = read_model_atmosphere("data/sun.mod")
    #    atm = Korg.PlanarAtmosphere(atm.layers[40:45])

    #    wls = Korg.Wavelengths(15000, 15500)

    #    A_X = format_A_X()
    #    @test synthesize(atm, [], A_X, 15000, 15500).wavelengths ≈ wls
    #    @test synthesize(atm, [], A_X, 15000, 15500; air_wavelengths=true).wavelengths ≈
    #          Korg.air_to_vacuum.(wls)
    #    @test_throws ArgumentError synthesize(atm, [], A_X, 15000, 15500; air_wavelengths=true,
    #                                          wavelength_conversion_warn_threshold=1e-20)
    #    @test_throws ArgumentError synthesize(atm, [], A_X, 2000, 8000, air_wavelengths=true)

    #    # test multiple line windows
    #    r1 = 5000:0.01:5001
    #    r2 = 6000:0.01:6001
    #    sol1 = synthesize(atm, [], A_X, r1; hydrogen_lines=true)
    #    sol2 = synthesize(atm, [], A_X, [r2]; hydrogen_lines=true)
    #    sol3 = synthesize(atm, [], A_X, [r1, r2]; hydrogen_lines=true)

    #    @test sol1.wavelengths == sol3.wavelengths[sol3.subspectra[1]]
    #    @test sol2.wavelengths == sol3.wavelengths[sol3.subspectra[2]]
    #    @test sol1.flux == sol3.flux[sol3.subspectra[1]]
    #    @test sol2.flux == sol3.flux[sol3.subspectra[2]]
    #end
end