# tests of things in src/utils.jl, not utilities for testing.  That's test/utilities.jl

@testset "utils" begin
    @testset "merge_bounds" begin
        mbounds = Korg.merge_bounds([(1, 3), (2, 4), (5, 6)], 0)
        @test mbounds[1] == [(1, 4), (5, 6)]
        @test mbounds[2] == [[1, 2], [3]]

        mbounds = Korg.merge_bounds([(1, 3), (2, 6), (5, 7)], 0)
        @test mbounds[1] == [(1, 7)]
        @test mbounds[2] == [[1, 2, 3]]

        mbounds = Korg.merge_bounds([(2, 6), (5, 7), (1, 3)], 0)
        @test mbounds[1] == [(1, 7)]
        @test mbounds[2] == [[3, 1, 2]]
    end

    @testset "LSF" begin
        wls = 5900:0.35:6100
        R = 1800.0
        spike_ind = 286 #in the middle
        flux = zeros(length(wls))
        flux[spike_ind] = 5.0

        # should't matter if you split wls into 2 ranges
        @test Korg.compute_LSF_matrix([5900:0.35:5935, 5935.35:0.35:6100], wls, R) ≈
              Korg.compute_LSF_matrix(wls, wls, R)

        # should't matter if R is a function
        @test Korg.compute_LSF_matrix(wls, wls, wl -> R) ≈ Korg.compute_LSF_matrix(wls, wls, R)

        # synth_wls can be a slightly fuzzy grid
        obs_wls = 5000:1.0:5010
        range_wls = 5000:0.01:5010
        array_wls = collect(range_wls) # range -> array
        array_wls[2:end-1] .+= 1e-8 * ones(length(range_wls) - 2) # slighly perturb (except the endpoints)
        @test Korg.compute_LSF_matrix(array_wls, obs_wls, R) ≈
              Korg.compute_LSF_matrix(range_wls, obs_wls, R)

        convF = Korg.apply_LSF(flux, wls, R)
        convF_5sigma = Korg.apply_LSF(flux, wls, R; window_size=5)
        convF_mat = Korg.compute_LSF_matrix(wls, wls, R) * flux
        convF_mat5 = Korg.compute_LSF_matrix(wls, wls, R; window_size=5) * flux
        convF_mat_vec = Korg.compute_LSF_matrix([5900:0.35:5935, 5935.35:0.35:6100], wls, R) * flux
        convF_changing_R = Korg.apply_LSF(flux, wls, wl -> wl / 6000 * R)
        convF_mat_changing_R = Korg.compute_LSF_matrix(wls, wls, wl -> wl / 6000 * R) * flux

        @test convF_mat_changing_R≈convF_changing_R rtol=1e-10

        downsampled_wls = 5950:0.4:6050
        convF_mat_downsample = Korg.compute_LSF_matrix(wls, downsampled_wls, R; window_size=5) *
                               flux

        # normalized?
        @test sum(flux)≈sum(convF) rtol=1e-3
        @test sum(flux)≈sum(convF_5sigma) rtol=1e-3
        @test sum(flux)≈sum(convF_mat) rtol=1e-3
        @test sum(flux)≈sum(convF_mat5) rtol=1e-3
        @test sum(flux)≈sum(convF_mat5) rtol=1e-3
        @test sum(flux)≈sum(convF_mat_vec) rtol=1e-3
        @test sum(flux)≈sum(convF_changing_R) rtol=1e-3
        @test sum(flux)≈sum(convF_mat_changing_R) rtol=1e-3
        @test sum(flux)*step(wls)≈sum(convF_mat_downsample)*step(downsampled_wls) rtol=1e-3

        # preserves line center?
        @test argmax(convF) == spike_ind
        @test argmax(convF_5sigma) == spike_ind
        @test argmax(convF_mat) == spike_ind
        @test argmax(convF_mat5) == spike_ind
        @test argmax(convF_mat_vec) == spike_ind
        @test argmax(convF_changing_R) == spike_ind
        @test argmax(convF_mat_changing_R) == spike_ind

        # make sure the default window_size values are OK
        @test assert_allclose(convF, convF_mat5; atol=1e-4)
        @test assert_allclose(convF_5sigma, convF_mat5; rtol=1e-11)
        @test assert_allclose(convF_mat, convF_mat5; atol=1e-4)

        # but also check that they are definitely doing something
        @test !(convF ≈ convF_5sigma)
        @test !(convF_mat ≈ convF_mat5)
    end

    @testset "rotation" begin
        # this implementation is less accurate than the one in Korg, but it produces correct-ish results
        function naive_apply_rotation(flux, wls::R, vsini, ε=0.6) where R<:AbstractRange
            vsini *= 1e5 # km/s to cm/s
            newFtype = promote_type(eltype(flux), eltype(wls), typeof(vsini), typeof(ε))
            newF = zeros(newFtype, length(flux))

            c1 = 2(1 - ε)
            c2 = π * ε / 2

            # step(wls) makes things normalized on the grid, and the factor of v_L in Gray becomes Δλrot
            # (because we are working in wavelenths) and moves inside the loop
            denominator = π * (1 - ε / 3) / step(wls)

            for i in 1:length(flux)
                Δλrot = wls[i] * vsini / Korg.c_cgs
                nwls = Int(floor(Δλrot / step(wls)))
                window = max(1, i - nwls):min(length(flux), i + nwls)

                x = (wls[i] .- wls[window]) ./ Δλrot
                one_less_x2 = @. 1 - x^2

                @. newF[window] .+= flux[i] * (c1 * sqrt(one_less_x2) + c2 * one_less_x2) /
                                    (denominator * Δλrot)
            end
            newF
        end

        wls = Korg.Wavelengths((4090, 5010))
        flux = zeros(length(wls))
        flux[990:1010] .= 1

        @testset for vsini in [0.0, 1e-10, 1.0, 5.0, 10.0, 20.0]
            @testset for ε in [0.1, 0.6, 0.9]
                # also test handling of multiple wl ranges
                @testset for wls in [wls, Korg.Wavelengths([4090:0.01:5007, 5007.01:0.01:5010])]
                    rflux = Korg.apply_rotation(flux, wls, vsini, ε)

                    # rotational kernel is normalized
                    @test sum(flux)≈sum(rflux) rtol=1e-2

                    if vsini > 1.0 && (length(wls.wl_ranges) == 1)
                        @test assert_allclose_grid(rflux,
                                                   naive_apply_rotation(flux, wls.wl_ranges[1],
                                                                        vsini,
                                                                        ε),
                                                   [("λ", wls * 1e8, "Å")]; atol=1e-2,
                                                   print_rachet_info=false)
                    end
                end
            end
        end
    end

    @testset "air <--> vacuum" begin
        wls = collect(2000.0:π:10000.0)
        @test Korg.vacuum_to_air.(Korg.air_to_vacuum.(wls))≈wls rtol=1e-3
        @test Korg.air_to_vacuum.(Korg.vacuum_to_air.(wls))≈wls rtol=1e-3

        #try it in cgs
        wls .*= 1e8
        @test Korg.vacuum_to_air.(Korg.air_to_vacuum.(wls))≈wls rtol=1e-3
        @test Korg.air_to_vacuum.(Korg.vacuum_to_air.(wls))≈wls rtol=1e-3

        #units should be automatically chosen
        @test Korg.vacuum_to_air.(Korg.air_to_vacuum.(wls * 1e-8) * 1e8)≈wls rtol=1e-3
        @test Korg.air_to_vacuum.(Korg.vacuum_to_air.(wls * 1e-8) * 1e8)≈wls rtol=1e-3
        @test Korg.vacuum_to_air.(Korg.air_to_vacuum.(wls) * 1e8)*1e-8≈wls rtol=1e-3
        @test Korg.air_to_vacuum.(Korg.vacuum_to_air.(wls) * 1e8)*1e-8≈wls rtol=1e-3
    end
end
