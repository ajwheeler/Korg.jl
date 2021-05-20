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
    elseif length(toks) == 1
        "_I"
    else
        throw(ArgumentError("invalid species code"))
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

#This type represents an individual line.
struct Line{F} 
    wl::F                     #given in nm, convert to cm
    log_gf::F                 #unitless
    species::String           
    E_lower::F                #eV
    gamma_rad::F              #s^-1
    gamma_stark::F            #s^-1
    vdW::Union{F, Tuple{F,F}} #either log(Γ_vdW) per electron or (σ, α) from ABO theory

    """
    Construct a `Line` with a possibly packed vdW parameter (sigma.alpha) format.  If vdW < 0,
    interpret it as log10(Γ) per particle.  Otherwise, interpret it as packed ABO parameters.
    """
    function Line(wl::F, log_gf::F, species::String, E_lower::F, gamma_rad::F, gamma_stark::F,
                  vdW::F) where F <: Real
        new{F}(wl, log_gf, species, E_lower, gamma_rad, gamma_stark, 
               if vdW > 0
                   floor(vdW) * bohr_radius_cgs * bohr_radius_cgs, vdW - floor(vdW)
               else 
                   10^vdW
               end
              )
    end
end
"""
Construct a `Line` without explicit broadening parameters.
"""
function Line(wl::F, log_gf::F, species::String, E_lower::F) where F <: Real
    Line(wl, log_gf, species, E_lower, approximate_radiative_gamma(wl, log_gf),
         approximate_gammas(wl, species, E_lower)...)
end

#pretty-print lines in REPL and jupyter notebooks
function Base.show(io::IO, m::MIME"text/plain", line::Line)
    print(io, line.species, " ", round(line.wl*1e8, digits=6), " Å")
end

"""
    approximate_radiative_gamma(wl, log_gf)

Approximate radiate broadening parameter.
"""
function approximate_radiative_gamma(wl, log_gf) 
    e = electron_charge_cgs
    m = electron_mass_cgs
    c = c_cgs
    8π^2 * e^2 / (m * c * wl^2) * 10^log_gf
end

"""
A simplified form of the Unsoeld (1995) approximation for van der Waals and Stark broadening at 
10,000 K.  The stark broadening 
Returns log10(γ_vdW).

In the calculation of n*², uses the approximation that
\\overbar{r}^2 = 5/2 {n^*}^4 / Z^2
which neglects the dependence on the angular momentum quantum number, l, in the the form given by
Warner 1967.

Used for atomic lines with no vdW and stark broadening info in the line list.
"""
function approximate_gammas(wl, species, E_lower; ionization_energies=ionization_energies)
    if ismolecule(species)
        return 0.0,0.0
    end

    ionization = split(species, '_')[2]
    Z = if ionization == "I"
        1
    elseif ionization == "II"
        2
    elseif ionization == "III"
        3
    end
    χ = ionization_energies[strip_ionization(species)][Z]
    c = c_cgs
    h = hplanck_eV
    k = kboltz_cgs
    E_upper = E_lower + (h * c / wl)

    nstar4_upper = (Z^2 * RydbergH_eV / (χ - E_upper))^2
    #From Cowley 1971
    γstark = 0.77e-18 * nstar4_upper * wl^2

    Δrbar2 = (5/2) * RydbergH_eV^2 * Z^2 * (1/(χ - E_upper)^2 - 1/(χ - E_lower)^2)
    if χ < E_upper
        println("Warning: for the $(species) line at $(Int(floor(wl*1e8))), the upper energy level"*
                " exceeds the ionization energy (E_upper) > $(χ)). Using null broadening params.")
        γvdW = 0.0
    else
        #From R J Rutten's course notes. An equivalent form can be found in Gray 2005.
        γvdW = 6.33 + 0.4log10(Δrbar2) + 0.3log10(10_000) + log10(k)
    end

    γstark, γvdW
end

"""
    new_line_imputing_zeros(wl, log_gf, species, E_lower, gamma_rad, gamma_stark, vdW)

Construct a new line treating broadening params equal to 0 as missing (how VALD represents missing
values).
"""
function new_line_imputing_zeros(wl, log_gf, species, E_lower, gamma_rad, gamma_stark, vdW)
    if gamma_rad == 0
        gamma_rad = approximate_radiative_gamma(wl, log_gf)
    end
    if (gamma_stark == 0) || (vdW == 0)
        approx_stark, approx_vdW = approximate_gammas(wl, species, E_lower)
        gamma_stark += (gamma_stark == 0)*approx_stark
        vdW += (vdW == 0)*approx_vdW
    end
    Line(wl, log_gf, species, E_lower, gamma_rad, gamma_stark, vdW)
end

"""
    read_line_list(fname; format="kurucz")

Parse the provided line list. in "Kurucz" format.
Pass `format="kurucz"` for a [Kurucz line list](http://kurucz.harvard.edu/linelists.html),
`format="vald"` for a Vald line list, and `format="moog"` for a MOOG line list.

Note that dissociation energies in a MOOG line list will be ignored.
"""
function read_line_list(fname::String; format="vald") :: Vector{Line}
    format = lowercase(format)
    linelist = open(fname) do f
        if format == "kurucz"
            parse_kurucz_linelist(f)
        elseif format == "vald"
            parse_vald_linelist(f)
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
        new_line_imputing_zeros(
            parse(Float64, line[1:11])*1e-7,
            parse(Float64, line[12:18]),
            parse_species_code(strip(line[19:24])),
            min(E_levels...),
            expOrZero(parse(Float64, line[81:86])),
            expOrZero(parse(Float64, line[87:92])),
            parse(Float64, line[93:98]))
    end
end

function _vald_to_korg_species_code(s)
     symbol, num = split(s[2:end-1], ' ')
     num = parse(Int, num)
     symbol * "_" * numerals[num]
end

function parse_vald_linelist(f)
    lines = collect(eachline(f))

    #figure out how big the header is
    firstline = findfirst(lines) do line
        line[1] == '\''
    end

    #air or vacuum wls?
    wl_header = split(lines[firstline-1])[3]
    if contains(wl_header, "air")
        wl_transform = air_to_vacuum
    elseif contains(wl_header, "vac")
        wl_transform = identity
    else
        throw(ArgumentError(
            "Can't parse line list.  I don't understant this wavelength column name: " * wl_header))
    end

    #vald short or long format?
    if isuppercase(lines[firstline][2]) && isuppercase(lines[firstline+1][2])
        Δ = 1 # short format
        shortformat = true
    else 
        Δ = 4 #long format
        shortformat = false
    end
    lines = lines[firstline:Δ:end]
    lastline = -1 + findfirst(lines) do line
        !((line[1] == '\'') && isuppercase(line[2]))
    end
    lines = lines[1:lastline]

    #filter out ions beyond III
    filter!(lines) do line
        findfirst(split(_vald_to_korg_species_code(split(line, ',')[1]), '_')[2] .== numerals) <= 3
    end

    map(lines) do line
        toks = split(line, ',')
        if shortformat
            new_line_imputing_zeros(
                parse(Float64, toks[2])*1e-8,
                parse(Float64, toks[5]),
                _vald_to_korg_species_code(toks[1]),
                parse(Float64, toks[3]),
                expOrZero(parse(Float64, toks[6])),
                expOrZero(parse(Float64, toks[7])),
                parse(Float64, toks[8]))
        else
            new_line_imputing_zeros(
                parse(Float64, toks[2])*1e-8,
                parse(Float64, toks[3]),
                _vald_to_korg_species_code(toks[1]),
                parse(Float64, toks[4]),
                expOrZero(parse(Float64, toks[11])),
                expOrZero(parse(Float64, toks[12])),
                parse(Float64, toks[13]))
        end
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
