module SSSynth

    _data_dir = joinpath(@__DIR__, "../data") 
    include("constants.jl")
    include("partition_func.jl")
    include("saha_boltzmann.jl")
    include("linelist.jl")

    #load data when the package is imported. We might as well do this until we ship with alternative 
    #datasets
    ionization_energies = setup_ionization_energies()
    partition_funcs = setup_atomic_partition_funcs()

end # module
