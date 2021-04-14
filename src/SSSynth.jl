module SSSynth
    export synthesize

    _data_dir = joinpath(@__DIR__, "../data") 

    include("constants.jl")      #physical constants
    include("atomic_data.jl")    #symbols and atomic weights
    include("line_opacity.jl")   #opacity, line profile, voigt function
    include("partition_func.jl") #approximate partition functions
    include("statmech.jl")       #statistical mechanics, molecular equilibrium
    include("linelist.jl")       #parse line lists
    include("atmosphere.jl")     #parse model atmospheres

    #load data when the package is imported. We might as well do this until we ship with alternative 
    #datasets
    ionization_energies = setup_ionization_energies()
    partition_funcs = setup_partition_funcs()
    equilibrium_constants = setup_equilibrium_constants()

    include("continuum_opacity/continuum_opacity.jl") #Define continuum opacity functions.
    include("synthesize.jl")                          #solve radiative transfer equation
    include("LSF.jl")

end # module
