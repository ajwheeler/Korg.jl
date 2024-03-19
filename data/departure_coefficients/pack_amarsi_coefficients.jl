using CSV, DataFrames, ProgressMeter, HDF5

function identify_file(files, keyword)
    fs = filter(files) do filename
        occursin(keyword, filename)
    end
    if length(fs) == 0
        error("No file found containing $keyword")
    elseif length(fs) > 1
        error("Multiple files found containing $keyword")
    end
    fs[1]
end

directory = ARGS[1]
outfile = ARGS[2]
println("packing departure coefficients from $directory, storing in $outfile")
files = readdir(directory)
label_file = joinpath(directory, identify_file(files, "label"))
grid_file = joinpath(directory, identify_file(files, "grid"))
atmos_file = joinpath(directory, identify_file(files, "atmos"))
println("labels: $label_file")
println("grid:   $grid_file")
println("atmos:  $atmos_file")

model_atom_levels = DataFrame(species=String[], config=String[], term=String[], J=Float64[], E=Float64[])
for line in eachline(label_file)
    if line[1] == '#'
        continue
    end
    vals = (strip(line[6:9]),
            strip(line[11:40]),
            strip(line[41:51]),
            parse(Float64,line[52:57]),
            parse(Float64, line[58:end]))
    push!(model_atom_levels, vals)
end

println("reading departure coefficients")
b_table = CSV.read(grid_file, DataFrame, header=["atm_ind", "depth_ind", "X_Fe", "level_ind", "b"],
                   skipto=2, delim=" ", ignorerepeated=true)

println("reading atmosphere info")
atm_table = CSV.read(atmos_file, DataFrame, header=["atm_ind", "Teff", "logg", "Fe_H"],
                     skipto=2, delim=" ", ignorerepeated=true)

vals = map([:Teff, :logg, :Fe_H]) do col
    sort(unique(atm_table[!, col]))
end
index_dicts = map(vals) do vs
    Dict(vs .=> 1:length(vs))
end
lengths = length.(index_dicts)

atm_indices = map(eachrow(atm_table)) do row
    row.atm_ind => getindex.(index_dicts, [row.Teff, row.logg, row.Fe_H])
end |> Dict

x_fe_vals = sort(unique(b_table.X_Fe))
nx_fe = length(x_fe_vals)
x_fe_indices = Dict(x_fe_vals .=> 1:nx_fe)

n_levels = maximum(b_table.level_ind) #TODO get elsewhere?

b_array = ones(Float32, 56, n_levels, lengths..., nx_fe)
println(sizeof(b_array) / 1024^3, " GB")
@showprogress desc="re-packing coefs into array" for row in eachrow(b_table)
    b_array[row.depth_ind, row.level_ind, atm_indices[row.atm_ind]..., x_fe_indices[row.X_Fe]] = row.b
end

h5open(outfile, "w") do file
    write(file, "b_array", b_array)
    write(file, "Teff", vals[1])
    write(file, "logg", vals[2])
    write(file, "m_H", vals[3])
    write(file, "x_Fe", x_fe_vals)
    write(file, "species", model_atom_levels.species)
    write(file, "configuration", model_atom_levels.config)
    write(file, "term", model_atom_levels.term)
    write(file, "J", model_atom_levels.J)
    write(file, "E", model_atom_levels.E)
end