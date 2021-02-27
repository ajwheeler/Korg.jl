module SSSynth

    _data_dir = joinpath(@__DIR__, "../data") 
    include("constants.jl")
    include("partition_func.jl")
    include("saha_boltzmann.jl")

end # module
