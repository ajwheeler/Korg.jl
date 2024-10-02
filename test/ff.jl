@testset "Peach ff absorption" begin
    # Tables IV-XXIV of Peach 1970 lists "absorption coefficient" from bound-free absorption and
    # the combined coefficient from both bound-free + free-free absorption per ion (for the .
    # This function returns an arbitrary subset from that table that is used to check the accuracy
    # of our free-free absorption calculation. Tabulated ref data from Peach (1970) is provided in 
    # the form (α includes stim emission):
    # k = α(ν, T, ndens_H_II, nₑ)/((1 - exp(-h*ν / (k*T)) * ndens_He_I)
    #
    # Specifically, this function returns data for He I ff absorption from Table IV. This returns
    # the absorption coefficient in cm² per He I atom (it's not a typo that this is per He I, not
    # per He II). The data excludes the correction for stimulated emission.
    #
    # There is way more data in the table that could be used for self-consistency checks, but the
    # format seems difficult to digitize, so we use only a few rows from the He II table.

    # the participating species and partition func, plus the wavelengths, temps, and k vals for 
    # He I ff from Peach. Testing heavier elements is more complicated because Peach's and our
    # partition funcs don't match, though in principal we could reconstruct hers by using the energy
    # levels in Table I of Peach 1970.
    He_I_ff_vals = let
        Ts = [1e4, 1.1e4, 1.2e4, 1.3e4, 1.4e4, 1.5e4]

        λ_cm = [7000.0, 3680.0] ./ 1e8

        bf_7000 = [1.442e-27, 1.689e-26, 1.314e-25, 7.469e-25, 3.315e-24, 1.207e-23]
        both_7000 = [1.665e-27, 2.019e-26, 1.630e-25, 9.616e-25, 4.432e-24, 1.676e-23]
        ff_7000 = both_7000 .- bf_7000
        bf_3680 = [2.878e-28, 3.337e-27, 2.577e-26, 1.455e-25, 6.421e-25, 2.327e-24]
        both_3680 = [3.221e-28, 3.845e-27, 3.061e-26, 1.784e-25, 8.129e-25, 3.043e-24]
        ff_3680 = both_3680 .- bf_3680
        out_ff = hcat(ff_7000, ff_3680)

        # at these temperatures, the partition function isn't defined. Let's define our own
        Us = Dict(Korg.species"He I" => T -> 1.0,
                  Korg.species"He II" => T -> 2.0,
                  Korg.species"He III" => T -> 1.0)

        Korg.species"He II", Us, λ_cm, Ts, out_ff
    end

    for (spec, Us, λ_cm, Ts, ref_k) in [He_I_ff_vals] #in case we want to test more later
        ν = Korg.c_cgs ./ λ_cm

        # compute α/(n(X+) * nₑ) by calculating α with both number densities set to 1.
        # here, X+ is X with one less electron, i.e. if X is He I, X+ is He II
        actual_α_div_nₑ_nHeII = zeros(length(Ts), length(ν))
        for (i, T) in enumerate(Ts)
            Korg.ContinuumAbsorption.positive_ion_ff_absorption!(view(actual_α_div_nₑ_nHeII, i, :),
                                                                 ν, T, Dict([spec => 1.0]), 1.0)
        end

        saha_RHS = map(Ts) do T
            # the ratio (n(X) * nₑ)/n(X+) is given by the Saha Equation.
            # Korg.saha_ion_weights(...)[1] gives n(X+)/n(X+). We set nₑ = 1.
            nₑ = 1.0
            atomic_num = Korg.get_atoms(spec)[1]
            wII, wIII = Korg.saha_ion_weights(T, nₑ, atomic_num, Korg.ionization_energies, Us)
            weights = [1.0, wII, wIII]
            weights[spec.charge+1] / weights[spec.charge]
        end
        stim_emission_factor = (1 .- exp.(-Korg.hplanck_cgs .* ν' ./ (Korg.kboltz_cgs .* Ts)))

        korg_k = actual_α_div_nₑ_nHeII .* saha_RHS ./ stim_emission_factor

        @test assert_allclose_grid(korg_k, ref_k, [("T", Ts, "K"), ("ν", ν, "Hz")];
                                   atol=0, rtol=0.005)
    end
end
