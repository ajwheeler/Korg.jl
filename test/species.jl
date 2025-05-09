@testset "species and formulae" begin
    @testset "species parsing" begin
        # Kurucz-style numerical codes
        @test Korg.species"01.00" == Korg.species"H"
        @test Korg.species"101.0" == Korg.species"H2 I"
        @test Korg.species"01.0000" == Korg.species"H I"
        @test Korg.species"02.01" == Korg.species"He II"
        @test Korg.species"02.1000" == Korg.species"He II"
        @test Korg.species"0608" == Korg.species"CO"
        @test Korg.species"0606" == Korg.species"C2 I"
        @test Korg.species"606" == Korg.species"C2 I"
        @test Korg.species"0608.00" == Korg.species"CO I"
        @test Korg.species"812.0" == Korg.species"MgO"
        @test Korg.species"822.0" == Korg.species"TiO"
        @test Korg.Species(1.0) == Korg.species"H"
        @test Korg.Species(101.0) == Korg.species"H2 I"
        @test Korg.Species(2.1) == Korg.species"He II"
        @test Korg.Species(2.01) == Korg.species"He II"

        @test Korg.species"010108" == Korg.species"10108" == Korg.species"H2O"
        @test Korg.species"060606" == Korg.species"60606" == Korg.species"C3"

        # The normal constructor (NOT string literal) must be used to test for failure
        # species which contain more than 2 tokens are invalid
        @test_throws ArgumentError Korg.Species("06.05.04")
        # Korg only goes up to uranium (Z=92)
        @test_throws Exception Korg.Species("93.01")

        #traditional-ish notation
        @test Korg.species"OOO" == Korg.species"O3"
        @test Korg.species"H 1" == Korg.species"H I"
        @test Korg.species"H     1" == Korg.species"H I"
        @test Korg.species"H_1" == Korg.species"H I"
        @test Korg.species"H.I" == Korg.species"H I"
        @test Korg.species"H I" == Korg.species"H I"
        @test Korg.species"H 2" == Korg.species"H II"
        @test Korg.species"H2" == Korg.species"HH I"
        @test Korg.species"H" == Korg.species"H I"
        @test Korg.species"C2H4" ==
              Korg.Species(Korg.Formula([0x01, 0x01, 0x01, 0x01, 0x06, 0x06]), 0)
        @test Korg.species"H+" == Korg.species"H II"
        @test Korg.species"OH+" == Korg.species"OH II"
        @test Korg.species"OH-" == Korg.Species(Korg.Formula("OH"), -1)

        # prevent constructing species with charges < -1.
        @test_throws ArgumentError Korg.Species("H -1")
        @test_throws ArgumentError Korg.Species("C2 -2")
    end

    @testset "couting atoms" begin
        for spec in Korg.Species.(["H2O", "H II", "C2-", "HHHH"])
            @test Korg.n_atoms(spec) == length(Korg.get_atoms(spec))
        end
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

    @testset "get_atom" begin
        @test Korg.get_atom(Korg.species"O I") == 0x08
        @test Korg.get_atom(Korg.species"O II") == 0x08
        @test_throws ArgumentError Korg.get_atom(Korg.species"CO")
    end
end
