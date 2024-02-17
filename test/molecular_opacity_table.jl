@testset "precomputed molecular opacities" begin
    # use 10 Å of the APOGEE linelist for testing
    linelist = Korg.get_APOGEE_DR17_linelist()
    filter!(linelist) do line
        15020 < line.wl * 1e8 < 15030
    end

    # split into water and non-water lines
    water_mask = [l.species == Korg.species"H2O" for l in linelist]
    linelist_less_water = linelist[.!water_mask]
    water_lines = linelist[water_mask]

    wls = 15020:0.01:15030

    table = Korg.precompute_molecular_cross_section(water_lines, wls; verbose=false)

    atm = interpolate_marcs(4000.0, 4.5)

    sol_without = synthesize(atm, linelist, format_A_X(), [wls]; verbose=false)
    sol_with = synthesize(atm, linelist_less_water, format_A_X(), [wls]; molecular_opacity_tables=[table], verbose=false)

    @test assert_allclose(sol_without.flux, sol_with.flux; rtol=1e-3)
end