@testset "Wavelengths" begin
    @testset "constructor" begin
        @test_throws "wavelengths must be non-empty" Korg.Wavelengths([])

        # Non-sorted and overlapping ranges should throw
        @test_throws ArgumentError Korg.Wavelengths([5000:1.0:5010, 5005:1.0:5020])
        @test_throws ArgumentError Korg.Wavelengths([5020:1.0:5030, 5000:1.0:5010])

        # Single-element vector input
        wls_single = Korg.Wavelengths([12345.6])
        wls_single_range = Korg.Wavelengths(12345.6, 12345.6)
        @test wls_single == wls_single_range
        @test length(wls_single) == 1
        @test wls_single[1] ≈ 12345.6 * 1e-8

        wls = Korg.Wavelengths(15000, 15500)
        @test wls == Korg.Wavelengths((15000, 15500))
        @test wls == Korg.Wavelengths([(15000, 15500)])
        @test wls == Korg.Wavelengths(15000:0.01:15500)
        @test wls == Korg.Wavelengths([15000:0.01:15500])
        @test wls == Korg.Wavelengths(collect(15000:0.01:15500))

        # should also work if you pass cm instead of Å
        # ideally these should be the same to the bit
        @test wls ≈ Korg.Wavelengths((15000e-8, 15500e-8))
        @test wls ≈ Korg.Wavelengths([(15000e-8, 15500e-8)])
        @test wls == Korg.Wavelengths((15000:0.01:15500) * 1e-8)
        @test wls == Korg.Wavelengths([(15000:0.01:15500) * 1e-8])
        @test wls == Korg.Wavelengths(collect(15000:0.01:15500) * 1e-8)

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
            @test searchsortedfirst(wls, 3 * conversion_factor) == 2

            @test searchsortedlast(wls, 8 * conversion_factor) == 6
            @test searchsortedlast(wls, 4 * conversion_factor) == 2
            @test searchsortedlast(wls, 11 * conversion_factor) == 9

            wls = Korg.Wavelengths(1:1.0:10)
            @test searchsortedfirst(wls, 5.4 * conversion_factor) == 6
            @test searchsortedlast(wls, 5.6 * conversion_factor) == 5

            wls = Korg.Wavelengths([3:5, 11:0.5:12.5, 16:20])
            @test searchsortedfirst(wls, 2 * conversion_factor) == 1
            @test searchsortedfirst(wls, 11.9 * conversion_factor) == 6
            @test searchsortedfirst(wls, 45 * conversion_factor) == 13
            @test searchsortedlast(wls, 2 * conversion_factor) == 0
            @test searchsortedlast(wls, 4 * conversion_factor) == 2
            @test searchsortedlast(wls, 11 * conversion_factor) == 4
            @test searchsortedlast(wls, 13.1 * conversion_factor) == 7
            @test searchsortedlast(wls, 55 * conversion_factor) == 12
        end
    end

    @testset "array interface" begin
        wls = Korg.Wavelengths(7000, 7002, 1.0)
        @test length(wls) == 3
        @test size(wls) == (3,)
        @test wls[1] ≈ 7000 * 1e-8
        @test wls[2] ≈ 7001 * 1e-8
        @test wls[3] ≈ 7002 * 1e-8
    end

    @testset "eachwindow and eachfreq" begin
        wls = Korg.Wavelengths([8000:1.0:8002, 8100:2.0:8104])
        windows = collect(Korg.eachwindow(wls))
        @test windows == [(8000e-8, 8002e-8), (8100e-8, 8104e-8)]
        freqs = Korg.eachfreq(wls)
        @test issorted(freqs)
        @test length(freqs) == length(wls)
    end
end
