@testset "Fit" begin
    @testset "parameter scaling" begin
        params = (Teff = 3200, logg=4.5, m_H=-2, vmic=3.2, vsini=10, O=-1)
        sparams = Korg.Fit.scale(params)
        uparams = Korg.Fit.unscale(sparams)
        @test all(isapprox.(values(uparams), values(params); rtol=1e-3))
    end

    @testset "fit param validation" begin
        
    end
end