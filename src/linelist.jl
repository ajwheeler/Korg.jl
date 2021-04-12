const numerals = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"]

"""
    parse_species_code(code)

Parse the "species code" as it is often specified in line lists and return a the "astronomy" 
notation. `01.00` → `"H_I"`, `02.01` → `"He_II"`, `0608` → `"CO"`, etc.  
"""
function parse_species_code(code::AbstractString)
    if length(code) == 4 && !contains(code, ".")
        atomic_symbols[parse(Int, code[1:2])] * atomic_symbols[parse(Int, code[3:4])]
    elseif 4 <= length(code) <= 5
        Z, charge = split(code, '.')
        atomic_symbols[parse(Int, Z)] * "_" * numerals[parse(Int, charge)+1]
    else
        throw(ArgumentError("Invalid species code: " * code))
    end
end

"get the chemical symbol for the element of the species"
strip_ionization(code::AbstractString)::String = split(code, '_')[1]

ismolecule(species::String) = sum(isuppercase(c) for c in species) > 1 || '2' in species

"""
Get the atoms that make up a diatomic molecule
"""
function get_atoms(molecule)
    if molecule[end] == '2'
        el1 = molecule[1:end-1]
        el2 = molecule[1:end-1]
    elseif isuppercase(molecule[2])
        el1, el2 = molecule[1:1], molecule[2:end]
    elseif isuppercase(molecule[3])
        el1, el2 = molecule[1:2], molecule[3:end]
    else
        throw(ArgumentError("This doesn't look like a diatomic molecule: $(molecule)"))
    end
    el1, el2        
end

"""
    read_line_list(fname; format="kurucz")

Parse the provided line list. in "Kurucz" format.
Pass `format="kurucz"` for a [Kurucz line list](http://kurucz.harvard.edu/linelists.html).
Pass `format="vald"` for a Vald line list. 
"""
function read_line_list(fname::String; format="kurucz", skiplines=2) :: Vector{NamedTuple}
    #at present, we preserve whatever order the linelist is in (e.g. groups by species, sorted by
    #wavelength) because the line opacity code makes no assumptions about linelist order, but we may 
    #want to normalize it in the future.
    linelist = if format == "kurucz"
        linelist = open(fname) do f
            map(eachline(f)) do line
                #kurucz provides wavenumbers for "level 1" and "level 2", which is which is 
                #determined by parity
                E_levels = map((line[25:36], line[53:64])) do s
                    #abs because Kurucz multiplies predicted values by -1
                    abs(parse(Float64,s)) * c_cgs * hplanck_eV
                end

                (wl=parse(Float64, line[1:11])*1e-7,             #given in nm, convert to cm
                 log_gf=parse(Float64, line[12:18]),
                 species=parse_species_code(strip(line[19:24])),
                 E_upper=max(E_levels...),                       #eV
                 E_lower=min(E_levels...),                       #eV
                 log_gamma_rad=parse(Float64, line[81:86]),      #s^-1
                 log_gamma_stark=parse(Float64, line[87:92]),    #s^-1
                 log_gamma_vdW=parse(Float64, line[93:98]))      #s^-1 (van der Waals)
            end
        end
    elseif format == "vald"
        lines = open(fname) do f
            lines = collect(eachline(f))
            #take only every fourth line in the file skipping the header and footer
            #the other three file lines per spectral line don't contain info we need.
            lines[skiplines+1:4:findfirst(l->startswith(l, "*"), lines)-1]
        end

        map(lines) do line
                toks = split(line, ',')
                (wl=parse(Float64, toks[2])*1e-8,
                 log_gf=parse(Float64, toks[3]),
                 species=begin 
                     symbol, num = split(toks[1][2:end-1], ' ')
                     num = parse(Int, num)
                     symbol * "_" * numerals[num]
                 end,
                 E_upper=parse(Float64, toks[6]), 
                 E_lower=parse(Float64, toks[4]),
                 log_gamma_rad=parse(Float64, toks[end-3]),
                 log_gamma_stark=parse(Float64, toks[end-2]),
                 log_gamma_vdW=parse(Float64, toks[end-1])
                )
            end
    else
        throw(ArgumentError("$(format) is not a supported line list format"))
    end

    mask = map(linelist) do line
        split(line.species, "_")[2] in ["I", "II", "III"]
    end
    if sum(.! mask) > 0
        @info "omitting $(sum(.! mask)) lines of high (> III) ionization states"
    end
    linelist[mask]
end
