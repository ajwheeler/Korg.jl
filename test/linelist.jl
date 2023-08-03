@testset "linelists" begin 
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

        fname = "kurucz_cn.txt"
        kurucz_ll = read_linelist("data/linelists/"*fname, format="kurucz")
        @test issorted(kurucz_ll, by=l->l.wl)
        @test length(kurucz_ll) == 10
        @test kurucz_ll[1].wl ≈ 2.9262621445487408e-5
        @test kurucz_ll[1].log_gf == -7.204
        @test kurucz_ll[1].species == Korg.species"CN"
        @test kurucz_ll[1].E_lower ≈ 1.1177309389190437
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
        @test linelist[3].gamma_stark ≈ 0.00019044182974029873
        @test linelist[3].vdW == linelist[1].vdW

        @test linelist[4].gamma_rad == linelist[1].gamma_rad
        @test linelist[4].gamma_stark == linelist[3].gamma_stark
        @test linelist[4].vdW == 9.953360714197118e-8

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
        for (Z, isotopes) in Korg.isotopic_abundances
            @test sum(values(isotopes)) ≈ 1
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

    @testset "turbospectrum linelists" begin
        ll = read_linelist("data/linelists/Turbospectrum/goodlist"; format="turbospectrum") 
        @test ll[1].species == ll[3].species
        @test ll[1].log_gf != ll[3].log_gf
        @test ll[1].E_lower == ll[3].E_lower
        @test ll[1].gamma_rad == ll[3].gamma_rad
        @test ll[1].gamma_stark == ll[3].gamma_stark
        @test ll[1].vdW ≈ ll[3].vdW

        @test ll[2].species == ll[3].species
        @test ll[2].log_gf ≈ ll[3].log_gf
        @test ll[2].E_lower == ll[3].E_lower
        @test ll[2].gamma_rad == ll[3].gamma_rad
        @test ll[2].gamma_stark == ll[3].gamma_stark
        @test ll[2].vdW ≈ ll[3].vdW

        vac_ll = read_linelist("data/linelists/Turbospectrum/goodlist"; format="turbospectrum_vac") 
        for (l_air, l_vac) in zip(ll, vac_ll)
            # l_vac.wl is "really" an air wavelength, but it wasn't converted because we told Korg 
            # to read it in as vacuum
            @test l_air.wl ≈ Korg.air_to_vacuum(l_vac.wl) rtol=1e-8
        end

        @test_throws ErrorException read_linelist("data/linelists/Turbospectrum/badlines"; format="turbospectrum")
        @test_throws ErrorException read_linelist("data/linelists/Turbospectrum/badvdw"; format="turbospectrum")
    end
end
