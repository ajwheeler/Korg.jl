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
    for i in 2:length(χs)
        #I think this should get optimized away?  I'm open to not doing this, but I'm hoping it 
        #makes formulas easier to read
        mₑ = electron_mass_cgs 
        k = kboltz_cgs
        k_eV = kboltz_eV
        h = hplanck_cgs
        weights[i] = (weights[i-1]/nₑ * (Us[i](T)/Us[i-1](T)) * (2π*mₑ*k*T/h^2)^1.5 * 
                      exp(-χs[i-1]/k_eV/T))
    end
    weights ./ sum(weights)
end
