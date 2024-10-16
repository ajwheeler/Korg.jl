@testset "Wavelengths" begin
    @testset "constructor" begin
        wls = Korg.Wavelengths(15000, 15500)
        @test wls == Korg.Wavelengths((15000, 15500))
        @test wls == Korg.Wavelengths([(15000, 15500)])
        @test wls == Korg.Wavelengths(15000:0.01:15500)
        @test wls == Korg.Wavelengths([15000:0.01:15500])
        @test wls == Korg.Wavelengths(collect(15000:0.01:15500))

        @test Korg.Wavelengths(15000, 15500, 1.0) == Korg.Wavelengths(15000, 15500, 1)
        @test Korg.Wavelengths(15000, 15500, 1) == Korg.Wavelengths(15000:1.0:15500)

        # test automatic air to vacuum conversion
        vac_wls = 15000:0.01:15500
        air_wls = Korg.air_to_vacuum.(15000:0.01:15500)
        for wls in [
            Korg.Wavelengths(15000:0.01:15500; air_wavelengths=true).wl_ranges[1] * 1e8,
            Korg.Wavelengths([15000:0.01:15500]; air_wavelengths=true).wl_ranges[1] * 1e8,
            Korg.Wavelengths(collect(15000:0.01:15500); air_wavelengths=true).wl_ranges[1] * 1e8,
            Korg.Wavelengths(15000, 15500; air_wavelengths=true).wl_ranges[1] * 1e8,
            Korg.Wavelengths((15000, 15500); air_wavelengths=true).wl_ranges[1] * 1e8
        ]
            @test assert_allclose_grid(wls, air_wls, [vac_wls]; atol=1e-4, print_rachet_info=false)
        end

        @test_throws ArgumentError Korg.Wavelengths(15000:0.01:15500; air_wavelengths=true,
                                                    wavelength_conversion_warn_threshold=1e-20)
    end

    @testset "subspectrum_indices" begin
        ranges = [5000:1.0:5010, 5020:1.0:5031, 5040:1.0:5060]
        @assert Korg.subspectrum_indices(Korg.Wavelengths(ranges)) == [1:11, 12:23, 24:44]
    end

    @testset "searchsortedfirst/last" begin

        # give the target wavelengths in both Å and cm
        @testset for conversion_factor in [1, 1e-8]
            wls = Korg.Wavelengths(0.5 .+ (2:10))
            @test searchsortedfirst(wls, 4.0 * conversion_factor) == 3
            @test searchsortedfirst(wls, 1.5 * conversion_factor) == 1
            @test searchsortedfirst(wls, 3.0 * conversion_factor) == 2
            @test searchsortedlast(wls, 8 * conversion_factor) == 6
            @test searchsortedlast(wls, 4 * conversion_factor) == 2
            @test searchsortedlast(wls, 11 * conversion_factor) == 9

            wls = Korg.Wavelengths(1:1.0:10)
            @test searchsortedfirst(wls, 5.4 * conversion_factor) == 6
            @test searchsortedlast(wls, 5.6 * conversion_factor) == 5

            wls = Korg.Wavelengths([3:5, 11:0.5:12.5, 16:20])
            @test searchsortedfirst(wls, 2) == 1
            @test searchsortedfirst(wls, 11.9) == 6
            @test searchsortedfirst(wls, 45) == 13
            @test searchsortedlast(wls, 2) == 0
            @test searchsortedlast(wls, 4) == 2
            @test searchsortedlast(wls, 11) == 4
            @test searchsortedlast(wls, 13.1) == 7
            @test searchsortedlast(wls, 55) == 12
        end
    end
end
