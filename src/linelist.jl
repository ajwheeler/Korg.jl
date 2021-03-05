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


"""
    setup_line_list(fname)

Parse the provided line list in "Kurucz" format.
"""
function read_line_list(fname::String)
    open(fname) do f
        map(eachline(f)) do line
            (wl=parse(Float64, line[1:11]), 
             log_gf=parse(Float64, line[12:18]),
             species=parse_species_code(strip(line[19:24])),
             E=parse(Float64, line[25:36]),                  #cm^-1
             log_gamma_rad=parse(Float64, line[81:86]),      #s^-1
             log_gamma_stark=parse(Float64, line[87:92]),    #s^-1
             log_gamma_vdW=parse(Float64, line[93:98]))      #s^-1 (van der Waals)
        end
    end
end


