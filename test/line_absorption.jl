# tests of things in src/line_absorption.jl and src/hydrogen_line_absorption.jl

@testset "line profiles" begin
    @testset "generic line profile" begin
        Δ = 0.01
        wls = (4750:Δ:5250) * 1e-8
        Δ *= 1e-8
        amplitude = 7.0
        for σ in [1e-7, 1e-8, 1e-9], γ in [3e-8, 3e-9, 3e-10]
            ϕ = Korg.line_profile.(5e-5, σ, γ, amplitude, wls)
            # the profile isn't perfectly monotonic because the approximation has "seams" at v=5
            # this allows for slight nonmonotonicity
            @test all(diff(ϕ[1:Int(ceil(end / 2))]) .> -1e-3 * maximum(ϕ))
            @test all(diff(ϕ[Int(ceil(end / 2)):end]) .< 1e-3 * maximum(ϕ))
            @test 0.98 < sum(ϕ .* Δ) / amplitude < 1
        end
    end

    @testset "hydrogen stark profiles" begin
        # This test data was generated with Korg.hydrogen_line_absorption shortly
        # after writing the function. This data is consistent with the results
        # produced by the Fortran code distributed with Stehle & Hutcheon 1999
        fname = "data/lyman_absorption.h5"
        αs_ref = h5read(fname, "profile")

        fid = h5open(fname)
        T = HDF5.read_attribute(fid["profile"], "T")
        ne = HDF5.read_attribute(fid["profile"], "ne")
        nH_I = HDF5.read_attribute(fid["profile"], "nH_I")
        wls = Korg.Wavelengths((HDF5.read_attribute(fid["profile"], "start_wl"):
                                HDF5.read_attribute(fid["profile"], "wl_step"):
                                HDF5.read_attribute(fid["profile"], "stop_wl")))
        close(fid)

        αs = zeros(length(wls))
        Korg.hydrogen_line_absorption!(αs, wls, 9000.0, ne, nH_I, 0.0,
                                       Korg.default_partition_funcs[Korg.species"H_I"](log(9000.0)),
                                       0.0, 15e-7; use_MHD=false)
        @test assert_allclose_grid(αs_ref, αs, [("λ", wls * 1e8, "Å")]; atol=5e-9)

        #make sure that H line absorption doesn't return NaNs on inputs where it used to
        wls = Korg.Wavelengths(3800:0.01:4200)
        αs = zeros(length(wls))
        Korg.hydrogen_line_absorption!(αs, wls, 9000.0, 1.1e16, 1, 0.0,
                                       Korg.default_partition_funcs[Korg.species"H_I"](log(9000.0)),
                                       0.0, 15e-7)
        @assert all(.!isnan.(αs))
    end
end
