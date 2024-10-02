# this script generates log_polyatomic_equilibrium_constants.h5 using polyatomic_partition_funcs.h5
# and atomization_energies.csv 
# This is included in case we want to precompute equilibrium constants in the future, but is
# currently unused.

# atomic parition functions (and species mass calculations) from Korg are used
using HDF5, CSV, Korg

# load polyatomic partition functions
partition_funcs = h5open("polyatomic_partition_funcs.h5") do f
    map(f) do group
        spec = Korg.Species(HDF5.name(group)[2:end]) # cut off leading '/'
        Ts, Us = read(group["temp"]), read(group["partition_function"])
        spec, (Ts, Us)
    end |> Dict
end

# equlibrium constants for polyatomics (done here because the partition functions are used)
atomization_Es = CSV.File("atomization_energies.csv")

h5open("log_polyatomic_equilibrium_constants.h5", "w") do f
    for (spec, D00) in zip(Korg.Species.(atomization_Es.spec), atomization_Es.energy)
        D00 *= 0.01036 # convert from kJ/mol to eV
        Ts, Us = partition_funcs[spec]

        log_pKs = map(zip(Ts, Us)) do (T, U)
            Zs = Korg.get_atoms(spec)
            Us_ratio = (prod([Korg.default_partition_funcs[Korg.Species(Korg.Formula(Z), 0)](log(T))
                              for Z in Zs]) / U)
            masses_ratio = prod([Korg.atomic_masses[Z] for Z in Zs]) / Korg.get_mass(spec)
            translational_U_factor = (2π * Korg.kboltz_cgs * T / Korg.hplanck_cgs^2)^1.5
            K = translational_U_factor^(length(Zs) - 1) * masses_ratio^1.5 * Us_ratio *
                exp(-D00 / (Korg.kboltz_eV * T))
            log10(K * (Korg.kboltz_cgs * T)^(length(Zs) - 1))
        end

        write(f, string(spec) * "/T", Ts)
        write(f, string(spec) * "/log_pK", log_pKs)
    end
end
