@testset "abundances" begin
    @test (format_A_X()
           == format_A_X(0)
           == format_A_X(0, 0)
           == format_A_X(Dict{String,Float64}())
           == format_A_X(Dict{Int,Float64}())
           == format_A_X(0, Dict(1 => 0.0); solar_relative=true)
           == format_A_X(0, 0, Dict(1 => 0.0); solar_relative=true)
           == format_A_X(0, Dict("H" => 0.0); solar_relative=true)
           == format_A_X(0, Dict(1 => 12.0); solar_relative=false)
           == format_A_X(0, Dict("H" => 12.0); solar_relative=false))

    # make sure silly H abundances are caught
    @test_throws ArgumentError format_A_X(0.0, Dict("H" => 0); solar_relative=false)
    @test_throws ArgumentError format_A_X(0.0, Dict(1 => 0); solar_relative=false)
    @test_throws ArgumentError format_A_X(0.0, Dict("H" => 12); solar_relative=true)
    @test_throws ArgumentError format_A_X(0.0, Dict(1 => 12); solar_relative=true)

    @testset "consitency of format_A_X and get_alpha_H/get_metals_H" begin
        @test Korg.get_alpha_H(format_A_X(0.1)) ≈ 0.1
        @test Korg.get_alpha_H(format_A_X(0.0, 0.1)) ≈ 0.1
        @test Korg.get_alpha_H(format_A_X(-0.2)) ≈ -0.2
        @test Korg.get_alpha_H(format_A_X(-2, -5)) ≈ -5

        @test Korg.get_metals_H(format_A_X(0.1)) ≈ 0.1
        @test Korg.get_metals_H(format_A_X(-0.2)) ≈ -0.2
        @test Korg.get_metals_H(format_A_X(0.1, 0.5)) ≈ 0.1
        @test Korg.get_metals_H(format_A_X(-0.2, 0.5)) ≈ -0.2

        # test consistency of alpha treatment by flipping each element to Inf and checking the result
        @testset for Z in 3:Korg.MAX_ATOMIC_NUMBER
            A_X = format_A_X()
            A_X[Z] = Inf
            if Z in Korg.default_alpha_elements
                @test Korg.get_metals_H(A_X) ≈ 0
                @test Korg.get_alpha_H(A_X) == Inf
            else
                @test Korg.get_metals_H(A_X) == Inf
                @test Korg.get_alpha_H(A_X) ≈ 0
            end
        end

        @test Korg.get_metals_H(Korg.grevesse_2007_solar_abundances;
                                solar_abundances=Korg.grevesse_2007_solar_abundances) ≈ 0
        @test Korg.get_alpha_H(Korg.grevesse_2007_solar_abundances;
                               solar_abundances=Korg.grevesse_2007_solar_abundances) ≈ 0
    end

    @test format_A_X(1.1) != format_A_X(1.1, 0)
    @test format_A_X(1.1)[50] == format_A_X(1.1, 0)[50] == format_A_X(-1, -2, Dict(50 => 1.1))[50]

    @testset for metallicity in [0.0, 0.5], abundances in [Dict(), Dict(:C => 1.1)],
                 solar_relative in [true, false]

        A_X = format_A_X(metallicity, abundances;
                         solar_relative=solar_relative)

        #correct absolute abundances?
        if "C" in keys(abundances)
            if solar_relative
                @test A_X[6] ≈ Korg.default_solar_abundances[6] + 1.1
            else
                @test A_X[6] ≈ 1.1
            end
        end
        @test A_X[7:end] ≈ Korg.default_solar_abundances[7:end] .+ metallicity
        @test A_X[1:2] == Korg.default_solar_abundances[1:2]
    end
end
