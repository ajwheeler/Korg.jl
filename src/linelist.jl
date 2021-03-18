"""
    parse_species_code(code)

Parse the "species code" as it is often specified in line lists and return a the "astronomy" 
notation. `01.00` → `"H_I"`, `02.01` → `"He_II"`, `0608.00` → `"CO"`, etc.  
"""
function parse_species_code(code::AbstractString)
    numerals = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"]
    if 4 <= length(code) <= 5
        Z, charge = split(code, '.')
        atomic_symbols[parse(Int, Z)] * "_" * numerals[parse(Int, charge)+1]
    elseif length(code) == 7
        throw(ArgumentError("Molecular species codes not yet supported"))
    else
        throw(ArgumentError("Invalid species code: " * code))
    end
end

#I should probably clean up my nomenclature here and throughout
"get the chemical symbol for the element of the species"
get_elem(code::AbstractString) = split(code, '_')[1]

"""
    read_line_list(fname)

Parse the provided line list in "Kurucz" format, 
[documented here](http://kurucz.harvard.edu/linelists.html).
"""
function read_line_list(fname::String) :: Vector{NamedTuple}
    linelist = open(fname) do f
        map(eachline(f)) do line
            (wl=parse(Float64, line[1:11])*10.,              #given in nm, convert to Å
             log_gf=parse(Float64, line[12:18]),
             species=parse_species_code(strip(line[19:24])),
             wavenumber=abs(parse(Float64, line[25:36])),    #cm^-1 (negative values are predicted)
             log_gamma_rad=parse(Float64, line[81:86]),      #s^-1
             log_gamma_stark=parse(Float64, line[87:92]),    #s^-1
             log_gamma_vdW=parse(Float64, line[93:98]))      #s^-1 (van der Waals)
        end
    end
    mask = map(linelist) do line
        split(line.species, "_")[2] in ["I", "II", "III"]
    end
    if sum(.! mask) > 0
        @warn "omitting $(sum(.! mask)) lines of high (> III) ionization states"
    end
    linelist[mask]
end


