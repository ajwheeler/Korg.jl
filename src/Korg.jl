module Korg
    export synthesize, read_linelist, read_model_atmosphere

    _data_dir = joinpath(@__DIR__, "../data") 

    include("constants.jl")                #physical constants
    include("atomic_data.jl")              #symbols and atomic weights
    include("linelist.jl")                 #parse linelists, define Line type
    include("line_opacity.jl")             #opacity, line profile, voigt function
    include("read_statmech_quantities.jl") #approximate Us, Ks, chis
    include("statmech.jl")                 #statistical mechanics, molecular equilibrium
    include("atmosphere.jl")               #parse model atmospheres
    include("transfer.jl")                 #radiative transfer integrals

    #load data when the package is imported. 
    const ionization_energies = setup_ionization_energies()
    const partition_funcs = setup_partition_funcs()
    const equilibrium_constants = setup_equilibrium_constants()
    const hline_stark_profiles = setup_hydrogen_stark_profiles()

    include("continuum_absorption/ContinuumAbsorption.jl")  #Define continuum absorption functions.
    include("synthesize.jl")                                #solve radiative transfer equation
    include("utils.jl")                                     #functions to apply LSF,
                                                            #vac<->air wls, etc

end # module
