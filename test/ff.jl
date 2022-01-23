using Korg
include("utilities.jl")

function _Peach70_TableIV_He_I_ff_ref()
    # Tables IV-XXIV of Peach 1970 lists "absorption coefficient" from bound-free absorption and
    # the combined coefficient from both bound-free + free-free absorption per ion (for the .
    # This function returns an arbitrary subset from that table that is used to check the accuracy
    # of our free-free absorption calculation
    #
    # Specifically, this function returns data for He I ff absorption from Table IV. This returns
    # the absorption coefficient in cm² per He I atom (it's not a typo that this is per He I, not
    # per He II). The data excludes the correction for stimulated emission.
    #
    # There is way more data in the table that could be used for self-consistency checks. But the
    # format seems difficult to digitize.

    # I should have selected higher temperature values (where the ff contribution is higher)
    T_row = [1e4, 1.1e4, 1.2e4, 1.3e4, 1.4e4, 1.5e4]
    
    bf_7000 = [ 1.442e-27, 1.689e-26, 1.314e-25, 7.469e-25, 3.315e-24, 1.207e-23]
    both_7000 = [1.665e-27, 2.019e-26, 1.630e-25, 9.616e-25, 4.432e-24, 1.676e-23]
    ff_7000 = both_7000 .- bf_7000


    bf_3680 = [2.878e-28, 3.337e-27, 2.577e-26, 1.455e-25, 6.421e-25, 2.327e-24]
    both_3680 = [3.221e-28, 3.845e-27, 3.061e-26, 1.784e-25, 8.129e-25, 3.043e-24]
    ff_3680 = both_3680 .- bf_3680

    out_ff = hcat(ff_7000, ff_3680)
    out_λ_cm = [7000.0, 3680.0] ./ 1e8
    out_λ_cm, T_row, out_ff
end

begin
    # tabulated ref data from Peach (1970) is provided in the form (α includes stim emission):
    # k = α(ν, T, ndens_H_II, nₑ)/((1 - exp(-h*ν / (k*T)) * ndens_He_I)
    λ_cm, T_arr, ref_k = _Peach70_TableIV_He_I_ff_ref()

    ν = Korg.c_cgs ./ λ_cm

    # let's compute α/(ndens_He_II * nₑ). We can do this by specifying ndens_He_II = 1 & nₑ = 1
    actual_α_div_nₑ_nHeII = Korg.ContinuumAbsorption._He_I_ff.(ν', T_arr, 1.0, 1.0)


    # the ratio (ndens_He_II * nₑ)/ndens_He_I is given by the Saha Equation
    # Korg.saha_ion_weights(...)[1] technically gives this divided by nₑ. We can get the desired
    # result by setting nₑ = 1.
    atomic_num = 2
    _nₑ = 1.0
    # at these temperatures, the partition function isn't defined. Let's define our own
    my_partition_funcs = Dict(Korg.Species("He I") => T->1.0,
                              Korg.Species("He II") => T->2.0,
                              Korg.Species("He III") => T->1.0)

    saha_RHS = map(T_arr) do cur_T
        Korg.saha_ion_weights(cur_T, _nₑ, atomic_num, Korg.ionization_energies,
                              my_partition_funcs)[1]
    end
    stim_emission_factor = (1 .- exp.(-Korg.hplanck_cgs .* ν' ./ (Korg.kboltz_cgs .* T_arr)))

    actual_k = actual_α_div_nₑ_nHeII .* saha_RHS ./ stim_emission_factor

    assert_allclose_grid(actual_k, ref_k, [("T", T_arr, "K"), ("ν", ν, "Hz")];
                         atol = 0, rtol = 0.01)
end
