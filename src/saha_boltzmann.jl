"""
   setup_ionization_energies([filename])

Parses the table of ionization energies and returns it as a dictionary mapping elements to
their ionization energies, `[χ₁, χ₂, χ₃]`.
"""
function setup_ionization_energies(fname=joinpath(_data_dir, 
                                                  "BarklemCollet2016-ionization_energies.dat"))
    open(fname, "r") do f
        d = Dict{String, Vector{Float64}}()
        for line in eachline(f)
            if line[1] != '#'        
                toks = split(strip(line))
                #the first token is the atomiz number, which we ignore
                d[toks[2]] = parse.(Float64, toks[3:end])
            end
        end
        d
    end
end

"""
    saha(χs, Us, T, nₑ)

Calculates the relative densities of the first n ionization states with ionization energies `χs` and 
partition functions `Us`, given the temperature, `T` [K], and electron density, `nₑ` [cm^-3].

Returns a 3-element vector of the same length which sums to 1.
"""
function saha(χs, Us, T, nₑ)
    weights = Vector{Float64}(undef, length(Us)) #there's probably a cleaner way to do this
    weights[1] = 1.
    for i in 2:length(Us)
        k_eV = kboltz_eV
        weights[i] = (2.0*weights[i-1]/nₑ * (Us[i](T)/Us[i-1](T)) *
                      free_partition_factor(electron_mass_cgs, T) * exp(-χs[i-1]/k_eV/T))
    end
    weights ./ sum(weights)
end

"""
The statitical multiplicities from the free movement of a particle. TODO reword.
Used in the Saha equation.

arguments
- `m` is the particle mass
- `T` is the temperature in K
"""
function free_partition_factor(m, T)
     k = kboltz_cgs
     h = hplanck_cgs
     (2π*m*k*T/h^2)^1.5
end

#using NLsolve

function molecular_equilibrium_equations(species)
    #remove ionization info
    species = get_elem.(species)

    #compute list of molecules, and the atoms of which they are composed
    molecules = collect(filter(ismolecule, species))
    atoms_pairs = map(collect(molecules)) do m
        #all molecules are diatomic, extract the two elements
        if isuppercase(m[2])
            el1, el2 = m[1:1], m[2:end]
        elseif isuppercase(m[3])
            el1, el2 = m[1:2], m[3:end]
        else
            @error "This doesn't look like a diatomic molecule: $(m)"
        end
        el1, el2        
    end
    atoms = collect(Set([first.(atoms_pairs) ; last.(atoms_pairs)]))

    #construct a mapping from each atom and molecular to a integer
    #these will be the indeces of each number density in the system equations, the xᵢs in the vector
    #for of the residual equation, F(x) = 0.
    var_indices = Dict{String, Int}()
    i = 1
    for a in atoms
        var_indices[i] = a
        i += 1
    end
    for m in molecules
        var_indeces[i] = m
    end

    #the residuals of the molecular equilibrium equations
    function residuals!(F, x)
        #molecular saha equations
    end
end


function molecular_equilibrium()


end
