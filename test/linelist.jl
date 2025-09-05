@testset "linelists" begin
    @testset "copy constructor" begin
        l1 = Korg.Line(5000.0, 0.0, Korg.species"Fe I", 1.0)
        @test l1 == Korg.Line(l1)
        @test Korg.Line(l1; wl=5001.0) != l1

        @test Korg.Line(l1; wl=5001.0).wl == 5001.0 * 1e-8
        @test Korg.Line(l1; log_gf=-0.5).log_gf == -0.5
        @test Korg.Line(l1; species=Korg.species"Fe II").species == Korg.species"Fe II"
        @test Korg.Line(l1; E_lower=2.0).E_lower == 2.0
        @test Korg.Line(l1; gamma_rad=1e-8).gamma_rad == 1e-8
        @test Korg.Line(l1; gamma_stark=1e-8).gamma_stark == 1e-8
        @test Korg.Line(l1; vdW=-8.0).vdW[1] == 1e-8
    end

    @testset "built-in linelists" begin
        # iterate over fns, not lists, because they make the output of the test suite way too long

        sol = nothing # cache chemical equilib
        @testset for linelist_fn in [Korg.get_VALD_solar_linelist,
            Korg.get_APOGEE_DR17_linelist,
            Korg.get_GALAH_DR3_linelist,
            Korg.get_GES_linelist
        ]
            linelist = linelist_fn()

            # test copy constructor
            @test linelist == Korg.Line.(linelist)

            @test issorted(linelist, by=l -> l.wl)

            # truncate model atmosphere for speed
            atm = Korg.read_model_atmosphere("data/sun.mod")
            atm = Korg.PlanarAtmosphere(atm.layers[1:3], atm.reference_wavelength)

            # make sure things run (types have caused problems in the past)
            λ = linelist[1].wl * 1e8
            synthesize(atm, linelist, format_A_X(), (λ, λ); use_chemical_equilibrium_from=sol)

            # test saving and loading
            filename = tempname() * ".h5"
            Korg.save_linelist(filename, linelist[1:1000])
            linelist2 = Korg.read_linelist(filename)
            @test linelist[1:1000] == linelist2
        end
    end

    @test_throws ArgumentError Korg.read_linelist("data/linelists/gfallvac08oct17.stub.dat";
                                                  format="abc")

    @testset "wls in either cm or Å" begin
        @test Korg.Line(5000.0, 0.0, Korg.species"Fe I", 1.0) ==
              Korg.Line(5e-5, 0.0, Korg.species"Fe I", 1.0)
    end

    @testset "ExoMol linelist parsing" begin
        # we pass this to load_ExoMol_linelist many times
        exomol_args = [
            Korg.species"CaH",
            "data/linelists/ExoMol/40Ca-1H__XAB_abridged.states",
            "data/linelists/ExoMol/40Ca-1H__XAB_abridged.trans"
        ]

        @testset "region specification" begin
            # whole thing
            linelist = Korg.load_ExoMol_linelist(exomol_args..., 6800, 6810; verbose=false)
            @test length(linelist) == 284
            @test all(6800 .<= [l.wl * 1e8 for l in linelist] .<= 6810)

            # no lines in this range
            linelist = Korg.load_ExoMol_linelist(exomol_args..., 5500, 6000; verbose=false)
            @test length(linelist) == 0

            # restrict to 5000-5025 Å
            linelist = Korg.load_ExoMol_linelist(exomol_args..., 6800, 6804; verbose=false)
            @test 0 < length(linelist) < 284 # some lines
            @test all(6800 .<= [l.wl * 1e8 for l in linelist] .<= 6804) # right wavelengths
        end

        @testset "empty linelist" begin
            # this use to crash when verbose=true
            linelist = Korg.load_ExoMol_linelist(exomol_args..., 3000, 4000; verbose=true)
            @test isempty(linelist)
        end

        @testset "isotopic corrections" begin
            # use line_strength_cutoff=-Inf to get all lines to make it easier to compare linelists

            # Test default behavior (most abundant isotopes)
            linelist_default = Korg.load_ExoMol_linelist(exomol_args..., 3000, 9000; verbose=false,
                                                         line_strength_cutoff=-Inf)
            # Test with explicit isotopes parameter (most abundant)
            isotopes = [(20, 40), (1, 1)]  # 40Ca, 1H
            linelist_explicit = Korg.load_ExoMol_linelist(exomol_args..., 3000, 9000; isotopes,
                                                          line_strength_cutoff=-Inf, verbose=false)

            # Should be identical since we're using the same isotopes
            @test linelist_default == linelist_explicit

            # Test with different isotopes (should have different log_gf values)
            # D has a spin degeneracy of 3 (H is 2) and a abundance of 1e-10 (by fiat)
            wrong_isotopes = [(20, 40), (1, 2)]
            linelist_rare = Korg.load_ExoMol_linelist(exomol_args..., 3000, 9000;
                                                      line_strength_cutoff=-Inf,
                                                      isotopes=wrong_isotopes, verbose=false)
            log_gf_deltas = map(linelist_default, linelist_rare) do l1, l2
                l2.log_gf - l1.log_gf
            end
            @test all(log_gf_deltas .≈ log10(2 / 3 * 1e-10))
        end

        @testset "ExoMol linelist line strength filtering" begin
            # Test with different line strength cutoffs
            linelist_loose = Korg.load_ExoMol_linelist(exomol_args..., 6800, 6804;
                                                       line_strength_cutoff=-20, verbose=false)

            linelist_strict = Korg.load_ExoMol_linelist(exomol_args..., 6800, 6804;
                                                        line_strength_cutoff=-10, verbose=false)

            # Loose cutoff should have more lines
            @test length(linelist_loose) >= length(linelist_strict)

            # All lines in strict should also be in loose
            strict_wavelengths = Set([l.wl for l in linelist_strict])
            loose_wavelengths = Set([l.wl for l in linelist_loose])
            @test strict_wavelengths ⊆ loose_wavelengths

            # Test with different temperature for line strength evaluation
            linelist_cold = Korg.load_ExoMol_linelist(exomol_args..., 6800, 6804;
                                                      T_line_strength=2000.0, verbose=false)

            linelist_hot = Korg.load_ExoMol_linelist(exomol_args..., 6800, 6804;
                                                     T_line_strength=5000.0, verbose=false)

            # Different temperatures may result in different numbers of lines
            # (this is more of a smoke test to ensure the parameter works)
            @test length(linelist_cold) >= 0
            @test length(linelist_hot) >= 0
        end

        @testset "ExoMol linelist error handling" begin
            # Test error handling for invalid isotopes
            @test_throws KeyError Korg.load_ExoMol_linelist(exomol_args..., 6800, 6804;
                                                            isotopes=[(20, 999), (1, 1)],
                                                            verbose=false)

            # Test error handling for invalid atomic numbers
            @test_throws KeyError Korg.load_ExoMol_linelist(exomol_args..., 6800, 6804;
                                                            isotopes=[(999, 40), (1, 1)],
                                                            verbose=false)
        end
    end

    @testset "kurucz linelist parsing" begin
        for fname in ["gfallvac08oct17.stub.dat", "gfallvac08oct17-missing-col.stub.dat"]
            kurucz_ll = Korg.read_linelist("data/linelists/" * fname; format="kurucz")
            @test issorted(kurucz_ll, by=l -> l.wl)
            @test length(kurucz_ll) == 987
            @test kurucz_ll[1].wl ≈ 0.0007234041763337705
            @test kurucz_ll[1].log_gf == -0.826
            @test kurucz_ll[1].species == Korg.species"Be_II"
            @test kurucz_ll[1].E_lower ≈ 17.360339371573698
            @test kurucz_ll[1].gamma_rad ≈ 8.511380382023759e7
            @test kurucz_ll[1].gamma_stark ≈ 0.003890451449942805
            @test kurucz_ll[1].vdW[1] ≈ 1.2302687708123812e-7
        end

        @testset "kurucz molecular " begin
            fname = "kurucz_cn.txt"
            @test_throws ArgumentError kurucz_ll=Korg.read_linelist("data/linelists/" * fname;
                                                                    format="kurucz")
            #@test issorted(kurucz_ll, by=l->l.wl)
            #@test length(kurucz_ll) == 10
            #@test kurucz_ll[1].wl ≈ 2.9262621445487408e-5
            #@test kurucz_ll[1].log_gf == -7.204
            #@test kurucz_ll[1].species == Korg.species"CN"
            #@test kurucz_ll[1].E_lower ≈ 1.1177309389190437
        end
    end

    @testset "vald short format, ABO, missing params" begin
        linelist = Korg.read_linelist("data/linelists/linelist.vald")
        @test length(linelist) == 6
        @test linelist[1].wl ≈ 3000.0414 * 1e-8
        @test linelist[1].log_gf == -2.957
        @test linelist[1].species == Korg.species"Fe_I"
        @test linelist[1].E_lower ≈ 3.3014
        @test linelist[1].gamma_rad ≈ 1.905460717963248e7
        @test linelist[1].gamma_stark ≈ 0.0001230268770812381
        @test linelist[1].vdW[1] ≈ 4.6773514128719815e-8

        #test imputation of missing broadening parameters
        @test linelist[2].gamma_rad ≈ 818252.5391161365
        @test linelist[2].gamma_stark == linelist[1].gamma_stark
        @test linelist[2].vdW == linelist[1].vdW

        @test linelist[3].gamma_rad == linelist[2].gamma_rad
        @test linelist[3].gamma_stark ≈ 0.00019044182974029873
        @test linelist[3].vdW == linelist[1].vdW

        @test linelist[4].gamma_rad == linelist[1].gamma_rad
        @test linelist[4].gamma_stark == linelist[3].gamma_stark
        @test linelist[4].vdW[1] == 9.953360714197118e-8

        @test linelist[5].gamma_rad == linelist[1].gamma_rad
        @test linelist[5].gamma_stark == linelist[1].gamma_stark
        @test linelist[5].vdW == linelist[4].vdW

        #ABO params
        @test linelist[6].vdW[1] ≈ 1.3917417470792187e-14
        @test linelist[6].vdW[2] ≈ 0.227
    end

    @testset "vald various formats" begin
        short_all = Korg.read_linelist("data/linelists/short-extract-all.vald")
        long_all_cm_air = Korg.read_linelist("data/linelists/long-extract-all-air-wavenumber.vald")
        long_all_cm_air_noquotes = Korg.read_linelist("data/linelists/long-extract-all-air-wavenumber-noquotes.vald")
        short_stellar = Korg.read_linelist("data/linelists/short-extract-stellar.vald")
        long_stellar = Korg.read_linelist("data/linelists/long-extract-stellar.vald")

        @test (length(short_all) == length(short_stellar) == length(long_all_cm_air) ==
               length(long_all_cm_air_noquotes) == length(long_stellar) == 2)

        #these should be identical since there was no unit conversion
        @test long_all_cm_air[1] == long_all_cm_air_noquotes[1]
        @test short_all[1] == short_stellar[1] == long_stellar[1]

        #when there was unit conversion, they should be approximately equal
        @test short_all[1].wl ≈ long_all_cm_air[1].wl
        @test short_all[1].log_gf == long_all_cm_air[1].log_gf
        @test short_all[1].species == long_all_cm_air[1].species
        @test short_all[1].E_lower≈long_all_cm_air[1].E_lower atol=1e-3
        @test short_all[1].gamma_rad == long_all_cm_air[1].gamma_rad
        @test short_all[1].gamma_stark == long_all_cm_air[1].gamma_stark
        @test short_all[1].vdW == long_all_cm_air[1].vdW
    end

    @testset "vald isotopic scaling" begin
        short_all = Korg.read_linelist("data/linelists/isotopic_scaling/short-all-unscaled.vald")
        long_all = Korg.read_linelist("data/linelists/isotopic_scaling/long-all-unscaled.vald")
        short_stellar = Korg.read_linelist("data/linelists/isotopic_scaling/short-stellar-unscaled.vald")
        long_stellar = Korg.read_linelist("data/linelists/isotopic_scaling/long-stellar-unscaled.vald")
        scaled = Korg.read_linelist("data/linelists/isotopic_scaling/scaled.vald")
        unscaled_lists = [short_all, long_all, short_stellar, long_stellar]

        for list in unscaled_lists
            for field in [:wl, :species, :E_lower, :gamma_rad, :gamma_stark, :vdW]
                @test getfield.(scaled, field) == getfield.(list, field)
            end
            for (line1, line2) in zip(scaled, list)
                @test line1.log_gf.≈line2.log_gf rtol=1e-2
            end
        end

        #make sure the default isotopic abundances sum to 1.
        for (Z, isotopes) in Korg.isotopic_abundances
            @test sum(values(isotopes)) ≈ 1
        end
    end

    @testset "Kurucz isotopic scaling (#463)" begin
        # this linelist contains the HFS components of one line for a few different isotopes
        my_iso = Korg.read_linelist("data/linelists/gfallvac08oct17_ba";
                                    format="kurucz") # use my isotopic abundances
        his_iso = Korg.read_linelist("data/linelists/gfallvac08oct17_ba"; format="kurucz",
                                     isotopic_abundances=nothing) # pull from linelist

        # the components should sum to the log(gf) of the line with no HFS or isotopic splitting
        # this one recovers the right log(gf)
        @test log10(sum(10^l.log_gf for l in my_iso))≈-0.03 atol=0.01
        # this one is slightly off because Kurucz doesn't unclude many digits
        @test log10(sum(10^l.log_gf for l in his_iso))≈-0.01 atol=0.01
    end

    moog_linelist = Korg.read_linelist("data/linelists/s5eqw_short.moog"; format="moog")
    moog_linelist_as_air = Korg.read_linelist("data/linelists/s5eqw_short.moog"; format="moog_air")
    @testset "moog linelist parsing" begin
        @test all(Korg.air_to_vacuum(l1.wl) .≈ l2.wl
                  for (l1, l2) in zip(moog_linelist, moog_linelist_as_air))

        @test issorted(moog_linelist, by=l -> l.wl)
        @test moog_linelist[1].wl ≈ 3729.807 * 1e-8
        @test moog_linelist[1].log_gf ≈ -0.280
        @test moog_linelist[1].species == Korg.species"Ti_I"
        @test moog_linelist[2].E_lower ≈ 3.265

        # test isotope parsing
        @test moog_linelist[4].species == Korg.species"MgH"
        @test moog_linelist[4].log_gf ≈ 0.52 + log10(Korg.isotopic_abundances[12][24])
        @test moog_linelist[5].species == Korg.species"C2"
        @test moog_linelist[5].log_gf ≈
              -0.082 + log10(Korg.isotopic_abundances[6][12]) +
              log10(Korg.isotopic_abundances[6][13])
        @test moog_linelist[6].species == Korg.species"Mn I"
        @test moog_linelist[6].log_gf ≈ -3.363 + log10(Korg.isotopic_abundances[25][55])
    end

    @testset "turbospectrum linelists" begin
        ll = Korg.read_linelist("data/linelists/Turbospectrum/goodlist"; format="turbospectrum")
        @test ll[1].species == ll[3].species
        @test ll[1].log_gf != ll[3].log_gf
        @test ll[1].E_lower == ll[3].E_lower
        @test ll[1].gamma_rad == ll[3].gamma_rad
        @test ll[1].gamma_stark == ll[3].gamma_stark

        @test ll[1].vdW[1] ≈ ll[2].vdW[1]
        @test ll[1].vdW[2] ≈ ll[2].vdW[2]

        @test ll[2].species == ll[3].species
        @test ll[2].log_gf ≈ ll[3].log_gf
        @test ll[2].E_lower == ll[3].E_lower
        @test ll[2].gamma_rad == ll[3].gamma_rad
        @test ll[2].gamma_stark == ll[3].gamma_stark

        @test ll[3].vdW[1] ≈ 8.89802482263476e-7

        vac_ll = Korg.read_linelist("data/linelists/Turbospectrum/goodlist";
                                    format="turbospectrum_vac")
        for (l_air, l_vac) in zip(ll, vac_ll)
            # l_vac.wl is "really" an air wavelength, but it wasn't converted because we told Korg
            # to read it in as vacuum
            @test l_air.wl≈Korg.air_to_vacuum(l_vac.wl) rtol=1e-8
        end

        @test_throws ErrorException Korg.read_linelist("data/linelists/Turbospectrum/badlines";
                                                       format="turbospectrum")
    end
end
