# this script scrapes enthalpies of formation at 0K from CCCBDB at NIST, then uses them to calculate 
# atomization energies for all the exomol polyatomic molecules in polyatomic_partition_funcs.h5, 
# which it writes to atomization_energies.csv

using Gumbo, HTTP, Cascadia, DataFrames, CSV, HDF5, Korg
using DataStructures: OrderedDict

# grab the CCCBDB table of enthalpies of formation at 0K
# modified from https://gist.github.com/scls19fr/9ea2fd021d5dd9a97271da317bff6533
function parsehtmltable(hdoc)
    qs = eachmatch(Selector("table"), hdoc.root)  # Return an Array{HTMLNode,1}
    tables = DataFrame[]
    for helm_table in qs
        column_names = String[]
        d_table = OrderedDict{String,Vector{String}}()
        for (i, row) in enumerate(eachmatch(Selector("tr"), helm_table))
            if (i == 1)
                for (j, colh) in enumerate(eachmatch(Selector("th"), row))
                    colh_text = strip(nodeText(colh))
                    while (colh_text in column_names)  # column header must be unique
                        colh_text = colh_text * "_2"
                    end
                    push!(column_names, colh_text)
                end
                if column_names == []
                    break
                end
            else
                if (i == 2)
                    for colname in column_names
                        d_table[colname] = Vector{String}()
                    end
                end
                for (j, col) in enumerate(eachmatch(Selector("td"), row))
                    col_text = strip(nodeText(col))
                    colname = column_names[j]
                    push!(d_table[colname], col_text)
                end
            end
        end
        df = DataFrame(d_table)
        push!(tables, df)
    end
    tables
end
cccbdb_html = parsehtml(String(HTTP.get("https://cccbdb.nist.gov/hf0k.asp").body));
Hfg_table, cccbdb_refs = parsehtmltable(cccbdb_html)[4:5]
Hfg_table.spec = map(Hfg_table.Species) do sp_str
    try
        Korg.Species(sp_str)
    catch e
        missing # parsing fails because the species contains deuterium, or has more than 6 atoms
    end
end
filter!(Hfg_table) do row
    !ismissing(row.spec)
end
Hfg = map(eachrow(Hfg_table)) do row
    row.spec, parse(Float64, row."Hfg 0K")
end |> Dict

# all the polyatomics for which partition functions are  available
polyatomics = h5open("polyatomic_partition_funcs.h5") do f
    map(f) do group
        Korg.Species(HDF5.name(group)[2:end]) # cut off leading '/'
    end
end

# calculate atomization energies by summing enthalpies of formation at 0K and write them to a file
atomization_energies = DataFrame(; spec=Korg.Species[], energy=Float64[])
for spec in polyatomics
    mol_enthalpy = if spec in keys(Hfg)
        Hfg[spec]
    else
        println(spec, " has no enthalpy")
        continue
    end

    atomic_enthalpy_sum = 0
    for Z in Korg.get_atoms(spec)
        atom_spec = Korg.Species(Korg.Formula(Z), 0)
        if atom_spec in keys(Hfg)
            atomic_enthalpy_sum += Hfg[atom_spec]
        else
            NaN
        end
    end
    if isnan(atomic_enthalpy_sum)
        println(spec, " has an atom with no enthalpy")
        continue
    end
    push!(atomization_energies, (spec, -mol_enthalpy + atomic_enthalpy_sum))
end
CSV.write("atomization_energies.csv", atomization_energies)
