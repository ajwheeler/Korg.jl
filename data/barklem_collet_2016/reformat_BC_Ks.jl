using Korg, FITSIO, DataFrames, HDF5

charged_mol_constituents = FITS("/Users/wheeler.883/Dropbox/Korg/data/barklem_collet_2016/J_A+A_588_A96_table1.dat.fits") do f
    DataFrame(; mol=Korg.Species.(read(f[2], "MolId")),
              el1=Korg.Species.(read(f[2], "Molc1")),
              el2=Korg.Species.(read(f[2], "Molc2")))
end
filter!(charged_mol_constituents) do row
    row.mol.charge == 1
end

# these are all on the same tempurature grid
bclogKs = Korg.read_Barklem_Collet_table("BarklemCollet2016-equilibrium_constants.dat";
                                         transform=x -> x + 1) #convert from log(mks) to log(cgs))
filter!(bclogKs) do (spec, Ks)
    0 <= spec.charge <= 1
end
bclogKs = map(collect(keys(bclogKs))) do spec
    itp = bclogKs[spec]
    spec => (itp.t, itp.u)
end |> Dict

# the temps for each species are the same at this stage
lnTs = bclogKs[Korg.species"H2"][1]
temps = exp.(lnTs)

# my partition funcs don't go as low in T as the B&C Ks
lnTmask = lnTs .>= 0

for row in eachrow(charged_mol_constituents)
    Korg.get_atoms(row.el1) == Korg.get_atoms(row.el2) && continue

    charged_el, neutral_el = if row.el1.charge > 0
        Korg.get_atoms(row.el1)[1], Korg.get_atoms(row.el2)[1]
    else
        Korg.get_atoms(row.el2)[1], Korg.get_atoms(row.el1)[1]
    end

    # only continue if the lower Z is not the charged atom
    charged_el <= neutral_el && continue

    logKs = bclogKs[row.mol][2]

    # A = lower χ, B = higher χ
    # get parition functions of the neutral and singly ionized for of each element
    UA = Korg.default_partition_funcs[Korg.Species(Korg.Formula(neutral_el), 0)].(lnTs[lnTmask])
    UAp = Korg.default_partition_funcs[Korg.Species(Korg.Formula(neutral_el), 1)].(lnTs[lnTmask])
    UB = Korg.default_partition_funcs[Korg.Species(Korg.Formula(charged_el), 0)].(lnTs[lnTmask])
    UBp = Korg.default_partition_funcs[Korg.Species(Korg.Formula(charged_el), 1)].(lnTs[lnTmask])

    χA = Korg.ionization_energies[neutral_el][1]
    χB = Korg.ionization_energies[charged_el][1]

    U_fac = @. log10(UBp) + log10(UA) - log10(UB) - log10(UAp)
    χ_fac = @. ((χA - χB) / (Korg.kboltz_eV * temps[lnTmask]))

    lTs = copy(lnTs)
    lTs[.!lnTmask] .= -Inf
    logKs[.!lnTmask] .= NaN
    logKs[lnTmask] .+= U_fac .+ χ_fac

    bclogKs[row.mol] = (lTs, logKs)
end

ps = collect(pairs(bclogKs))
mols = string.(first.(ps))
lnTs = reduce(vcat, transpose.(first.(last.(ps))))
logKs = reduce(vcat, transpose.(last.(last.(ps))))

rm("barklem_collet_ks.h5") # doesn't error if the file doesn't exist
h5write("barklem_collet_ks.h5", "mols", mols)
h5write("barklem_collet_ks.h5", "lnTs", lnTs)
h5write("barklem_collet_ks.h5", "logKs", logKs)
