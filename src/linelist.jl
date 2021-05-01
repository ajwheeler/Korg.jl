const numerals = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"]

"""
    parse_species_code(code)

Parse the "species code" as it is often specified in line lists and return a the "astronomy" 
notation. 01.00 → "H_I", 02.01 → "He_II", 02.1000 → "He_II", 0608 → "CO_I", etc.  
"""
function parse_species_code(code::AbstractString)
    toks = split(code, '.')

    ionization = if length(toks) == 2
        if all(c == '0' for c in toks[2])
            "_I"
        else
            "_" * numerals[parse(Int, replace(toks[2], "0"=>""))+1]
        end
    else
        "_I"
    end

    atom_or_molecule = if length(toks[1]) < 3
        atomic_symbols[parse(Int, toks[1])]
    elseif length(toks[1]) == 4
        atomic_symbols[parse(Int, toks[1][1:2])] * atomic_symbols[parse(Int, toks[1][3:4])]
    end

    atom_or_molecule * ionization
end

"get the chemical symbol for the element of the species"
strip_ionization(code::AbstractString)::String = split(code, '_')[1]

"""
true if the string passed represents a molecule (with or without its ionization state)
"""
function ismolecule(species)
    count = 0
    for c in species
        if c == '_'
            break
        else
            count += isdigit(c) + isuppercase(c)
        end
    end
    count > 1
end

"""
Get the atoms that make up a diatomic molecule
"""
function get_atoms(molecule)
    if '_' in molecule
        molecule = split(molecule,'_')[1]
    end
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
This type represents an individual line.
"""
struct Line{F} 
    wl::F                     #given in nm, convert to cm
    log_gf::F                 #unitless
    species::String           
    E_lower::F                #eV
    gamma_rad::F              #s^-1
    gamma_stark::F            #s^-1
    vdW::Union{F, Tuple{F,F}} #either Γ_vdW per electron or (σ, α) from ABO theory

    """
    Construct a `Line` with a possibly packed vdW parameter (sigma.alpha) format.  If vdW < 0,
    interpret it as log10(Γ) per particle.  Otherwise, interpret it as packed ABO parameters.
    """
    function Line(wl::F, log_gf::F, species::String, E_lower::F, gamma_rad::F, gamma_stark::F,
                  vdW::F) where F <: AbstractFloat
        new{F}(wl, log_gf, species, E_lower, gamma_rad, gamma_stark, 
             if vdW > 0
                 (floor(vdW) * bohr_radius_cgs * bohr_radius_cgs, vdW - floor(vdW))
             elseif vdW == 0
                 0.0
             else
                 10^vdW
             end) 
    end
end
"""
Construct a `Line` without explicit broadening parameters.
"""
function Line(wl::F, log_gf::F, species::String, E_lower::F) where F <: AbstractFloat
    e = electron_charge_cgs
    m = electron_mass_cgs
    c = c_cgs
    Line(wl, log_gf, species, E_lower, 8π^2 * e^2 / (m * c * wl^2) * 10^log_gf, 0.0, 0.0)
end

"""
The Unsoeld (1995) approximation for van der Waals broadening.  Used for atomic lines with no vdW 
broadening info in the line list.
"""
function unsoeld(line; ionization_energies=ionization_energies)
    ionization = split(line.species, '_')[2]
    Z = if ionization == "I"
        1
    elseif ionization == "II"
        2
    elseif ionization == "III"
        3
    end
    χ = ionization_energies[strip_ionization(line)][Z]
    nstar2 = RydbergH_eV * Z*Z / (χ - line.E_lower)
    #rbar2 = 

    #10^()
end

#pretty-print lines in REPL and jupyter notebooks
function Base.show(io::IO, m::MIME"text/plain", line::Line)
    print(io, line.species, " ", round(line.wl*1e8, digits=6), " Å")
end

"""
    read_line_list(fname; format="kurucz")

Parse the provided line list. in "Kurucz" format.
Pass `format="kurucz"` for a [Kurucz line list](http://kurucz.harvard.edu/linelists.html),
`format="vald"` for a Vald line list, and `format="moog"` for a MOOG line list.

Note that dissociation energies in a MOOG line list will be ignored.
"""
function read_line_list(fname::String; format="kurucz", skiplines=2) :: Vector{Line}
    format = lowercase(format)
    linelist = open(fname) do f
        if format == "kurucz"
            parse_kurucz_linelist(f)
        elseif format == "vald"
            parse_vald_linelist(f, skiplines)
        elseif format == "moog"
            parse_moog_linelist(f)
        else
            throw(ArgumentError("$(format) is not a supported line list format"))
        end
    end

    mask = map(linelist) do line
        split(line.species, "_")[2] in ["I", "II", "III"]
    end
    if sum(.! mask) > 0
        @info "omitting $(sum(.! mask)) lines of high (> III) ionization states"
    end

    #ensure linelist is sorted
    if issorted(linelist[mask], by=l->l.wl)
        linelist[mask]
    else
        linelist = linelist[mask]
        sort!(linelist, by=l->l.wl)
        linelist
    end
end

#used in to parse vald and kurucz lineslists
expOrZero(x) = x == 0.0 ? 0.0 : 10.0^x

function parse_kurucz_linelist(f)
    map(eachline(f)) do line
        #kurucz provides wavenumbers for "level 1" and "level 2", which is which is 
        #determined by parity
        E_levels = map((line[25:36], line[53:64])) do s
            #abs because Kurucz multiplies predicted values by -1
            abs(parse(Float64,s)) * c_cgs * hplanck_eV
        end
        Line(parse(Float64, line[1:11])*1e-7,
             parse(Float64, line[12:18]),
             parse_species_code(strip(line[19:24])),
             min(E_levels...),
             expOrZero(parse(Float64, line[81:86])),
             expOrZero(parse(Float64, line[87:92])),
             parse(Float64, line[93:98]))
    end
end

function parse_vald_linelist(f, skiplines)
    lines = collect(eachline(f))
    #take only every fourth line in the file skipping the header and footer
    #the other three file lines per spectral line don't contain info we need.
    lines = lines[skiplines+1:4:findfirst(l->startswith(l, "*"), lines)-1]
    map(lines) do line
        toks = split(line, ',')
        Line(parse(Float64, toks[2])*1e-8,
             parse(Float64, toks[3]),
             begin 
                 symbol, num = split(toks[1][2:end-1], ' ')
                 num = parse(Int, num)
                 symbol * "_" * numerals[num]
             end,
             parse(Float64, toks[4]),
             expOrZero(parse(Float64, toks[11])),
             expOrZero(parse(Float64, toks[12])),
             parse(Float64, toks[13]))
    end
end

#todo support moog linelists with broadening parameters?
function parse_moog_linelist(f)
    lines = collect(eachline(f))
    #moog format requires blank first line
    linelist = map(lines[2:end]) do line
        toks = split(line)
        wl =
        parse(Float64, toks[4])
        Line(parse(Float64, toks[1]) * 1e-8, #convert Å to cm
             parse(Float64, toks[4]),
             parse_species_code(toks[2]),
             parse(Float64, toks[3]))
    end
    #TODO issue warning, don't autoconvert
    linelist
end

