@testset "Fit" begin
    @testset "parameter scaling" begin
        params = (Teff = 3200, logg=4.5, m_H=-2, vmic=3.2, vsini=10, O=-1)
        sparams = Korg.Fit.scale(params)
        uparams = Korg.Fit.unscale(sparams)
        @test all(isapprox.(values(uparams), values(params); rtol=1e-3))
    end

    @testset "fit param validation" begin
        @test_throws ArgumentError Korg.Fit.validate_params((; Teff = 3200), (;))
        @test_throws ArgumentError Korg.Fit.validate_params((; logg = 3200), (;))
        @test_throws ArgumentError Korg.Fit.validate_params((Teff=4500, logg = 3200, m_H = 0.1), (; m_H=0.1))

        _, fixed_params = Korg.Fit.validate_params((Teff=4500, logg=4.5), (;))
        @test fixed_params.m_H == 0
        @test fixed_params.vmic == 1
    end

    @testset "merge bounds" begin
        @test Korg.Fit.merge_bounds([(1, 3), (2, 4), (5,6)], 0) == [(1, 4), (5,6)]
        @test Korg.Fit.merge_bounds([(1, 3), (2, 6), (5,7)], 1.0) == [(1, 7)]
    end
end