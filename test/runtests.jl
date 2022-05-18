using Korg, Test, HDF5
import Dierckx

@testset "cubic spline interpolation" begin
    xs = 0:0.1:5
    ys = sin.(xs)

    ditp = Dierckx.Spline1D(xs, ys; k=3, bc="error", s=0.0)
    kitp = Korg.CubicSplines.CubicSpline(xs, ys)

    #these are only approximately the same because Dierckx is a BSpline library.
    x = 0:0.01:5
    @test Dierckx.evaluate(ditp, x) ≈ kitp.(x) rtol=1e-3

    @test_throws ArgumentError kitp(-1)
    @test_throws ArgumentError kitp(6)
end

@testset "atomic data" begin 
    @test (Korg.Natoms == length(Korg.atomic_masses) == length(Korg.asplund_2009_solar_abundances) 
            == length(Korg.asplund_2020_solar_abundances))
    @test (Korg.get_mass(Korg.Formula("CO")) ≈ 
           Korg.get_mass(Korg.Formula("C")) + Korg.get_mass(Korg.Formula("O")))
    @test Korg.get_mass(Korg.Formula("C2")) ≈ 2Korg.get_mass(Korg.Formula("C"))
end

@testset "ionization energies" begin
    @test length(Korg.ionization_energies) == 92
    @test Korg.ionization_energies[Korg.atomic_numbers["H"]] == [13.5984, -1.000, -1.000]
    @test Korg.ionization_energies[Korg.atomic_numbers["Ru"]] == [7.3605, 16.760, 28.470]
    @test Korg.ionization_energies[Korg.atomic_numbers["U"]] == [6.1940, 11.590, 19.800]
end

function _test_contained_slice(vals::AbstractVector, interval::Korg.Interval)

    idx = Korg.contained_slice(vals, interval)
    first_ind, last_ind = first(idx), last(idx)

    @assert first_ind >= 1 && last_ind <= length(vals)

    result = Korg.contained.(vals, Ref(interval))

    if all(result)
        @test (first_ind == 1) && (last_ind == length(vals))
    elseif any(result)
        @test last_ind >= first_ind
        @test all(.!result[1:first_ind-1])
        @test all(result[first_ind:last_ind])
        @test all(.!result[last_ind+1:length(vals)])
    else
        @test first_ind == last_ind+1
    end

end

@testset "Interval" begin

    # first make sure that the following cases are caught by the constructor:
    @test_throws AssertionError Korg.Interval(5,5)
    @test_throws AssertionError Korg.Interval(3,2)
    @test_throws AssertionError Korg.Interval(Inf,Inf)
    @test_throws AssertionError Korg.Interval(-Inf,-Inf)

    # check contained
    sample = Korg.Interval(3,10)
    @test !Korg.contained(3, sample)
    @test !Korg.contained(10, sample)
    @test Korg.contained(5.0, sample)
    @test Korg.contained(nextfloat(3.0), sample)
    @test Korg.contained(prevfloat(10.0), sample)

    # check contained_slice
    @testset "contained_slice" begin
        # we consider cases where the slice contains just a single element or multiple elements

        # first, try cases where everything is in-bounds
        _test_contained_slice([6.0], sample) # slice of 1 element
        _test_contained_slice([4.0, 5.0, 6.0, 7.0, 8.0, 9.0], sample) # slice of multiple elements

        # next, try cases where some values are out-of bounds
        _test_contained_slice([1.0, 6.0], sample) # slice of 1 element
        _test_contained_slice([1.0, 6.0, 12.0], sample)
        _test_contained_slice([6.0, 12.0], sample)

        _test_contained_slice([1.0, 2.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], sample)
        _test_contained_slice([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.5, 12.0], sample)
        _test_contained_slice([1.0, 2.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.5, 12.0], sample)

        # lastly, consider cases where all values are out of bounds
        _test_contained_slice([1.0], sample)
        _test_contained_slice([100.0], sample)

        _test_contained_slice([1.0, 2.0, 2.5], sample)
        _test_contained_slice([10.5, 12.5, 100.0], sample)
    end
end

include("continuum_absorption.jl") # test this after the "Interval" testset

"""
Compute nₑ (number density of free electrons) in a pure Hydrogen atmosphere, where `nH_tot` is the
total number density of H I and H II (in cm⁻³), the temperature is `T`, and `HI_partition_val` is
the the value of the H I partition function.

This is a relatively naive implementation. More numerically stable solutions exist.
"""
function electron_ndens_Hplasma(nH_tot, T, H_I_partition_val = 2.0)
    # Define the Saha equation as: nₑ*n_{H II} / n_{H I} = RHS
    # coef ∼ 4.829e15
    coef = 2.0 * (2.0*π*Korg.electron_mass_cgs*Korg.kboltz_cgs / Korg.hplanck_cgs^2)^1.5
    RHS = coef * T^1.5 * exp(-Korg.RydbergH_eV/(Korg.kboltz_eV*T))/H_I_partition_val
    # In a pure Hydrogen atmosphere: nₑ = n_{H II}. The Saha eqn becomes:  nₑ²/(nH_tot - ne) = RHS
    # We recast the Saha eqn as: a*nₑ² + b*nₑ + c = 0 and compute the coefficients
    a, b, c = (1.0, RHS, -1*RHS*nH_tot)
    # solve quadratic equation. Since b is always positive and c is always negative:
    #    (-b + sqrt(b²-4*a*c))/(2*a) is always ≥ 0
    #    (-b - sqrt(b²-4*a*c))/(2*a) is always negative
    nₑ = (-b + sqrt(b*b-4*a*c))/(2*a)
    nₑ
end

@testset "O I-III and CN partition functions are (nearly) monotonic in T" begin
    lnTs = 0:0.1:log(10_000.0)
    nearly_monotonic(Us) = all(diff(Us) .> -1e-4)
    @test nearly_monotonic(Korg.partition_funcs[Korg.species"O_I"].(lnTs))
    @test nearly_monotonic(Korg.partition_funcs[Korg.species"O_II"].(lnTs))
    @test nearly_monotonic(Korg.partition_funcs[Korg.species"O_III"].(lnTs))
    @test nearly_monotonic(Korg.partition_funcs[Korg.species"CN_I"].(lnTs))
end

@testset "stat mech" begin
    @testset "pure Hydrogen atmosphere" begin
        nH_tot = 1e15
        # specify χs and Us to decouple this testset from other parts of the code
        χs = [(Korg.RydbergH_eV, -1.0, -1.0)]
        Us = Dict([Korg.species"H_I"=>(T -> 2.0), Korg.species"H_II"=>(T -> 1.0)])
        # iterate from less than 1% ionized to more than 99% ionized
        for T in [3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 1e4, 1.1e4, 1.2e4, 1.3e4, 1.4e4, 1.5e5]
            nₑ = electron_ndens_Hplasma(nH_tot, T, 2.0)
            wII, wIII = Korg.saha_ion_weights(T, nₑ, 1, χs, Us)
            @test wIII ≈ 0.0 rtol = 1e-15
            rtol = (T == 1.5e5) ? 1e-9 : 1e-14
            @test wII/(1 + wII + wIII) ≈ (nₑ/nH_tot) rtol= rtol
        end
    end

    @testset "monotonic N ions Temperature dependence" begin
        weights = [Korg.saha_ion_weights(T, 1.0, 7, Korg.ionization_energies, 
                                            Korg.partition_funcs) for T in 1:100:10000]
        #N II + NIII grows with T === N I shrinks with T
        @test issorted(first.(weights) + last.(weights))
        
        # NIII grows with T
        @test issorted(last.(weights))
    end

    @testset "molecular equilibrium" begin
        #solar abundances
        abundances = Korg.get_absolute_abundances(0.0, Dict(), Korg.asplund_2020_solar_abundances, true)
        nₜ = 1e15 
        nₑ = 1e-3 * nₜ #arbitrary

        MEQs = Korg.molecular_equilibrium_equations(abundances, Korg.ionization_energies, 
                                                       Korg.partition_funcs, 
                                                       Korg.equilibrium_constants)
        @test MEQs.atoms == 0x01:Korg.Natoms

        n = Korg.molecular_equilibrium(MEQs, 5700.0, nₜ, nₑ)
        #make sure number densities are sensible
        @test (n[Korg.species"C_III"] < n[Korg.species"C_II"] < n[Korg.species"C_I"] < 
               n[Korg.species"H_II"] < n[Korg.species"H_I"])

        #total number of carbons is correct
        total_C = map(collect(keys(n))) do species
            if species == Korg.species"C2"
                n[species] * 2
            elseif 0x06 in Korg.get_atoms(species.formula)
                n[species]
            else
                0.0
            end
        end |> sum
        @test total_C ≈ abundances[Korg.atomic_numbers["C"]] * nₜ
    end
end

@testset "lines" begin
    @testset "linelists" begin 
        @testset "species codes" begin
            @test Korg.species"01.00"   == Korg.species"H"
            @test Korg.species"101.0"   == Korg.species"H2 I"
            @test Korg.species"01.0000" == Korg.species"H I"
            @test Korg.species"02.01"   == Korg.species"He II"
            @test Korg.species"02.1000" == Korg.species"He II"
            @test Korg.species"0608"    == Korg.species"CO"
            @test Korg.species"0606"    == Korg.species"C2 I"
            @test Korg.species"606"     == Korg.species"C2 I"
            @test Korg.species"0608.00" == Korg.species"CO I"
            @test Korg.species"OOO"     == Korg.species"O3"

            @test_throws ArgumentError Korg.Species("06.05.04")
            @test_throws Exception Korg.Species("99.01")
        end

        @testset "distinguish atoms from molecules" begin
            @test Korg.ismolecule(Korg.Formula("H2"))
            @test Korg.ismolecule(Korg.Formula("CO"))
            @test !Korg.ismolecule(Korg.Formula("H"))
            @test !Korg.ismolecule(Korg.Formula("Li"))
        end

        @testset "break molecules into atoms" begin
            @test Korg.get_atoms(Korg.Formula("CO")) == [0x06, 0x08]
            @test Korg.get_atoms(Korg.Formula("C2")) == [0x06, 0x06]
            @test Korg.get_atoms(Korg.Formula("MgO")) == [0x08, 0x0c]
        end

        @test_throws ArgumentError read_linelist("data/linelists/gfallvac08oct17.stub.dat";
                                                          format="abc")

        @testset "kurucz linelist parsing" begin
            for fname in ["gfallvac08oct17.stub.dat", "gfallvac08oct17-missing-col.stub.dat"]
                kurucz_ll = read_linelist("data/linelists/"*fname, format="kurucz")
                @test issorted(kurucz_ll, by=l->l.wl)
                @test length(kurucz_ll) == 987
                @test kurucz_ll[1].wl ≈ 0.0007234041763337705
                @test kurucz_ll[1].log_gf == -0.826
                @test kurucz_ll[1].species == Korg.species"Be_II"
                @test kurucz_ll[1].E_lower ≈ 17.360339371573698
                @test kurucz_ll[1].gamma_rad ≈ 8.511380382023759e7
                @test kurucz_ll[1].gamma_stark ≈ 0.003890451449942805
                @test kurucz_ll[1].vdW ≈ 1.2302687708123812e-7
            end
        end

        @testset "vald short format, ABO, missing params" begin
            linelist = read_linelist("data/linelists/linelist.vald")
            @test length(linelist) == 6
            @test linelist[1].wl ≈ 3000.0414 * 1e-8
            @test linelist[1].log_gf == -2.957
            @test linelist[1].species == Korg.species"Fe_I"
            @test linelist[1].E_lower ≈ 3.3014
            @test linelist[1].gamma_rad ≈ 1.905460717963248e7
            @test linelist[1].gamma_stark ≈ 0.0001230268770812381
            @test linelist[1].vdW ≈ 4.6773514128719815e-8

            #test imputation of missing broadening parameters
            @test linelist[2].gamma_rad ≈ 818252.5391161365
            @test linelist[2].gamma_stark == linelist[1].gamma_stark
            @test linelist[2].vdW == linelist[1].vdW

            @test linelist[3].gamma_rad == linelist[2].gamma_rad
            @test linelist[3].gamma_stark == 0.0001903497656604677
            @test linelist[3].vdW == linelist[1].vdW

            @test linelist[4].gamma_rad == linelist[1].gamma_rad
            @test linelist[4].gamma_stark == linelist[3].gamma_stark
            #this is failing on windows (due to https://github.com/JuliaLang/julia/issues/43088 ?)
            @test_skip linelist[4].vdW == 9.953360714197118e-8

            @test linelist[5].gamma_rad == linelist[1].gamma_rad
            @test linelist[5].gamma_stark == linelist[1].gamma_stark
            @test linelist[5].vdW == linelist[4].vdW

            #ABO params
            @test linelist[6].vdW[1] ≈ 1.3917417470792187e-14
            @test linelist[6].vdW[2] ≈ 0.227
        end

        @testset "vald various formats" begin
            short_all = read_linelist("data/linelists/short-extract-all.vald")
            long_all_cm_air = read_linelist("data/linelists/long-extract-all-air-wavenumber.vald")
            long_all_cm_air_noquotes = read_linelist("data/linelists/long-extract-all-air-wavenumber-noquotes.vald")
            short_stellar = read_linelist("data/linelists/short-extract-stellar.vald")
            long_stellar = read_linelist("data/linelists/long-extract-stellar.vald")

            @test (length(short_all) == length(short_stellar) == length(long_all_cm_air) == 
                   length(long_all_cm_air_noquotes) == length(long_stellar) == 2)

            #these should be identical since there was no unit conversion
            @test long_all_cm_air[1] == long_all_cm_air_noquotes[1]
            @test short_all[1] == short_stellar[1] == long_stellar[1]

            #when there was unit conversion, they should be approximately equal
            @test short_all[1].wl ≈ long_all_cm_air[1].wl
            @test short_all[1].log_gf == long_all_cm_air[1].log_gf
            @test short_all[1].species == long_all_cm_air[1].species
            @test short_all[1].E_lower ≈ long_all_cm_air[1].E_lower     atol=1e-3
            @test short_all[1].gamma_rad == long_all_cm_air[1].gamma_rad
            @test short_all[1].gamma_stark == long_all_cm_air[1].gamma_stark
            @test short_all[1].vdW == long_all_cm_air[1].vdW
        end

        @testset "vald isotopic scaling" begin
            short_all = read_linelist("data/linelists/isotopic_scaling/short-all-unscaled.vald")
            long_all = read_linelist("data/linelists/isotopic_scaling/long-all-unscaled.vald")
            short_stellar = read_linelist("data/linelists/isotopic_scaling/short-stellar-unscaled.vald")
            long_stellar = read_linelist("data/linelists/isotopic_scaling/long-stellar-unscaled.vald")
            scaled = read_linelist("data/linelists/isotopic_scaling/scaled.vald")
            unscaled_lists = [short_all, long_all, short_stellar, long_stellar]

            for list in unscaled_lists
                for field in [:wl, :species, :E_lower, :gamma_rad, :gamma_stark, :vdW]
                    @test getfield.(scaled, field) == getfield.(list, field)
                end
                for (line1, line2) in zip(scaled, list)
                    @test line1.log_gf .≈ line2.log_gf rtol=1e-2
                end
            end

            #make sure the default isotopic abundances sum to 1.
            isotopes = keys(Korg.isotopic_abundances)
            for Z in sort(collect(Set(first.(isotopes))))
                isos = filter(==(Z) ∘ first, isotopes)                                               
                @assert sum([Korg.isotopic_abundances[iso] for iso in isos]) ≈ 1                       
            end
        end

        moog_linelist = read_linelist("data/linelists/s5eqw_short.moog"; format="moog")
        @testset "moog linelist parsing" begin
            @test issorted(moog_linelist, by=l->l.wl)
            @test moog_linelist[1].wl ≈ 3729.807 * 1e-8
            @test moog_linelist[1].log_gf ≈ -0.280
            @test moog_linelist[1].species == Korg.species"Ti_I"
            @test moog_linelist[2].E_lower ≈ 3.265
        end
    end

    @testset "move_bounds" begin
        a = 0.5 .+ (1:9)
        for lb in [1, 3, 9], ub in [1, 5, 9]
            @test Korg.move_bounds(a, lb, ub, 5., 2.) == (3, 6)
            @test Korg.move_bounds(a, lb, ub, 0., 3.) == (1, 2)
            @test Korg.move_bounds(a, lb, ub, 6., 4.) == (2, 9)
            @test Korg.move_bounds(collect(a), lb, ub, 5., 2.) == (3, 6)
            @test Korg.move_bounds(collect(a), lb, ub, 0., 3.) == (1, 2)
            @test Korg.move_bounds(collect(a), lb, ub, 6., 4.) == (2, 9)
        end
    end

    @testset "line profile" begin
        Δ = 0.01
        wls = (4955 : Δ : 5045) * 1e-8
        Δ *= 1e-8
        amplitude = 7.0
        for Δλ_D in [1e-7, 1e-8, 1e-9], Δλ_L in [1e-8, 1e-9]
            ϕ = Korg.line_profile.(5e-5, 1/Δλ_D, Δλ_L, amplitude, wls)
            @test issorted(ϕ[1 : Int(ceil(end/2))])
            @test issorted(ϕ[Int(ceil(end/2)) : end], rev=true)
            @test 0.99 < sum(ϕ .* Δ)/amplitude < 1
        end
    end

    @testset "hydrogen stark profiles" begin
        # This test data was generated with Korg.hydrogen_line_absorption shortly
        # after writing the function. This data is consistent with the results
        # produced by the Fortran code distributed with Stehle & Hutcheon 1999
        fname = "data/lyman_absorption.h5"
        αs_ref = h5read(fname,  "profile")

        fid = h5open("data/lyman_absorption.h5") 
        T = HDF5.read_attribute(fid["profile"], "T")
        ne = HDF5.read_attribute(fid["profile"], "ne")
        nH_I = HDF5.read_attribute(fid["profile"], "nH_I")
        wls = (HDF5.read_attribute(fid["profile"], "start_wl") :
               HDF5.read_attribute(fid["profile"], "wl_step") : 
               HDF5.read_attribute(fid["profile"], "stop_wl") )
        close(fid)

        αs = zeros(length(wls))
        Korg.hydrogen_line_absorption!(αs, wls, 9000.0, 1e11, 1e13, 
                                       Korg.partition_funcs[Korg.species"H_I"](log(9000.0)), 
                                       Korg.hline_stark_profiles, 0.0)
        @test αs_ref ≈ αs rtol=1e-5
    end
end

@testset "atmosphere" begin
    @testset "plane-parallel atmosphere" begin
        #the MARCS solar model atmosphere
        atm = Korg.read_model_atmosphere("data/sun.mod")
        @test atm isa Korg.PlanarAtmosphere
        @test length(atm.layers) == 56
        @test issorted([l.temp for l in atm.layers])
        @test atm.layers[1].tau_5000 ≈ 0.00001209483645
        @test atm.layers[1].z == 6.931E+07
        @test atm.layers[1].temp == 4066.8
        @test atm.layers[1].electron_number_density ≈ 3.769664452210607e10
        @test atm.layers[1].number_density ≈ 4.75509171357701e14
    end
    @testset "spherical atmosphere" begin
        atm = Korg.read_model_atmosphere(
                "data/s6000_g+1.0_m0.5_t05_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
        @test atm isa Korg.ShellAtmosphere
        @test length(atm.layers) == 53 #drop layers hotter than 10^4 K
        @test issorted([l.temp for l in atm.layers])
        @test atm.R == 2.5827E+12
        @test atm.layers[1].tau_5000 ≈ 4.584584692493259e-5
        @test atm.layers[1].z == 2.222e11
        @test atm.layers[1].temp == 3935.2
        @test atm.layers[1].electron_number_density ≈ 1.7336231777439526e8
        @test atm.layers[1].number_density ≈ 1.5411190391302566e12
    end

    @testset "atmosphere type conversion" begin
        atm = Korg.read_model_atmosphere("data/sun.mod")
        atm2 = Korg.PlanarAtmosphere(Korg.ShellAtmosphere(atm, 7e10)) #arbitrary radius
        @test atm.layers == atm2.layers

        atm = Korg.read_model_atmosphere(
                "data/s6000_g+1.0_m0.5_t05_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
        atm2 = Korg.ShellAtmosphere(Korg.PlanarAtmosphere(atm), 1.0)
        @test [l.tau_5000 for l in atm.layers]                == [l.tau_5000 for l in atm2.layers]
        @test [l.z for l in atm.layers]                       == [l.z for l in atm2.layers]
        @test [l.temp for l in atm.layers]                    == [l.temp for l in atm2.layers]
        @test [l.number_density for l in atm.layers]          == [l.number_density for l in atm2.layers]
        @test [l.electron_number_density for l in atm.layers] == [l.electron_number_density for l in atm2.layers]
    end
end

@testset "synthesis" begin

    @testset "calculate absolute abundances" begin
        @test_throws ArgumentError Korg.get_absolute_abundances(0.0, Dict("H"=>13), 
                                                        Korg.asplund_2020_solar_abundances, true)

        @testset for metallicity in [0.0, 0.5], abundances in [Dict(), Dict("C"=>1.1)], solar_relative in [true, false]
            nx_nt = Korg.get_absolute_abundances(metallicity, abundances, 
                                                Korg.asplund_2020_solar_abundances, solar_relative)
            #correct absolute abundances?
            if "C" in keys(abundances)
                if solar_relative
                    @test log10(nx_nt[6]/nx_nt[1]) + 12 ≈ Korg.asplund_2020_solar_abundances[6] + 1.1
                else
                    @test log10(nx_nt[6]/nx_nt[1]) + 12 ≈ 1.1
                end
            end
            @test log10(nx_nt[2]/nx_nt[1]) + 12 ≈ Korg.asplund_2020_solar_abundances[2]
            @test log10(nx_nt[Korg.atomic_numbers["Ba"]]/nx_nt[1]) + 12 ≈ 
                Korg.asplund_2020_solar_abundances[Korg.atomic_numbers["Ba"]] + metallicity

            #normalized?
            @test sum(nx_nt) ≈ 1
        end

    end

    @testset "trapezoid rule" begin
        #gaussian PDF should integral to 1.
        pdf(x) = exp(-1/2 * x^2) / sqrt(2π)
        xs = -10:0.1:10
        fs = pdf.(xs)

        partial_ints = similar(xs)
        Korg.cumulative_trapezoid_rule!(partial_ints, xs, fs)

        naive_partial_integrals = map(1:length(xs)) do i
            Korg.trapezoid_rule(xs[1:i], fs[1:i])
        end

        #do they match the analytic solution?
        @test Korg.trapezoid_rule(xs, fs) ≈ 1.0  atol=1e-5
        @test partial_ints[end] ≈ 1.0 atol=1e-5

        #do they match each other?
        @test naive_partial_integrals ≈ partial_ints
    end
end

@testset "LSF" begin
    wls = 5000:0.35:6000
    R = 1800.0
    flux = zeros(Float64, length(wls))
    flux[500] = 5.0

    convF = Korg.constant_R_LSF(flux, wls, R)
    #normalized?
    @test sum(flux) ≈ sum(convF)

    #preserves line center?
    @test argmax(convF) == 500
end

@testset "air <--> vacuum" begin
    wls = collect(2000.0:π:10000.0)
    @test Korg.vacuum_to_air.(Korg.air_to_vacuum.(wls)) ≈ wls rtol=1e-3
    @test Korg.air_to_vacuum.(Korg.vacuum_to_air.(wls)) ≈ wls rtol=1e-3

    #try it in cgs
    wls .*= 1e8
    @test Korg.vacuum_to_air.(Korg.air_to_vacuum.(wls)) ≈ wls rtol=1e-3
    @test Korg.air_to_vacuum.(Korg.vacuum_to_air.(wls)) ≈ wls rtol=1e-3

    #units should be automatically chosen
    @test Korg.vacuum_to_air.(Korg.air_to_vacuum.(wls*1e-8)*1e8) ≈ wls rtol=1e-3
    @test Korg.air_to_vacuum.(Korg.vacuum_to_air.(wls*1e-8)*1e8) ≈ wls rtol=1e-3
    @test Korg.vacuum_to_air.(Korg.air_to_vacuum.(wls)*1e8)*1e-8 ≈ wls rtol=1e-3
    @test Korg.air_to_vacuum.(Korg.vacuum_to_air.(wls)*1e8)*1e-8 ≈ wls rtol=1e-3
end

using SpecialFunctions: expint
@testset "transfer" begin
    xs = 2:0.01:8
    @test Korg.exponential_integral_2.(xs) ≈ expint.(2, xs) rtol=1e-3
end

@testset "synthesize wavelength handling" begin
    atm = read_model_atmosphere("data/sun.mod")
    wls = 15000:0.01:15500
    @test synthesize(atm, [], 15000, 15500).wavelengths ≈ wls
    @test synthesize(atm, [], 15000, 15500; air_wavelengths=true).wavelengths ≈ Korg.air_to_vacuum.(wls)
    @test_throws ArgumentError synthesize(atm, [], 15000, 15500; air_wavelengths=true, 
                                          wavelength_conversion_warn_threshold=1e-20)
    @test_throws ArgumentError synthesize(atm, [], 2000, 8000, air_wavelengths=true)
end

@testset "autodiff" begin
    using ForwardDiff

    for atm_file in ["data/sun.mod",
             "data/s6000_g+1.0_m0.5_t05_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod"]
        atm = read_model_atmosphere(atm_file)
        linelist = read_linelist("data/linelists/5000-5005.vald")
        wls = 6564:0.01:6565
        flux(p) = synthesize(atm, linelist, wls; metallicity=p[1], abundances=Dict(["Ni"=>p[2]]), 
                             vmic=p[3]).flux
        #make sure this works.
        ∇f = ForwardDiff.jacobian(flux, [0.0, 0.0, 1.5])
    end
end
