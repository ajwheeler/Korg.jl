import Base: (==), hash

"""
Represents an atom or molecule, irespective of its charge.
"""
struct Baryon
    atoms::Vector{String}

    """
        Baryon(code::String)

        Construct a Baryon from an encoded string form.  This can be a MOOG-style numeric code, i.e.
        "0801" for OH, or an atomic or molecular symbol, i.e. "FeH", "Li", or "C2".
    """
    function Baryon(code::AbstractString)
        #handle numeric codes, e.g. 0801 -> OH
        if all(isdigit(c) for c in code)
            if length(code) <= 2
                return Baryon(parse(Int,code))
            elseif length(code) <= 4
                if length(code) == 3  
                    code = "0"*code
                end
                return Baryon(parse(Int, code[1:2]), parse(Int, code[3:4]))
            else
                throw(ArgumentError("numeric codes for molecules with more than 4 chars are not "*
                                    "supported"))
            end
        end
        #otherwise, code should be "OH", "FeH", "Li", "C2", etc.
        inds = filter(1:length(code)) do i
            isdigit(code[i]) || isuppercase(code[i])
        end
        push!(inds, length(code)+1)
        subcode = map(1:(length(inds)-1)) do j
            code[inds[j]:inds[j+1]-1]
        end

        atoms = String[]
        for s in subcode
            num = tryparse(Int, s)
            if num isa Int
                previous = atoms[end]
                for _ in 1:(num-1)
                    push!(atoms, previous)
                end
            else
                push!(atoms, s)
            end
        end

        sort!(atoms, by=s->atomic_numbers[s])
        new(String.(atoms))
    end

    """
        Baryon(symbols::Int...)

    Create a Baryon by providing the atomic symbols of its components as arguments
    """
    function Baryon(symbols::AbstractString...)
        new(sort(collect(String.(symbols)), by=s->atomic_numbers[s]))
    end

    """
        Baryon(Z::Int...)

    Create a Baryon by providing the atomic numbers of its components as arguments
    """
    function Baryon(Z::Int...)
        new([atomic_symbols[z] for z in sort(collect(Z))])
    end

    """
        Baryon(symbols::Vector{AbstractString})

    A baryon composed of the atoms in `symbols`, which must be sorted by atomic number.
    """
    function Baryon(symbols::Vector{AbstractString})
        @assert issorted(symbols, by=s->atomic_numbers[s])
        new(String.(symbols))
    end
end

function (==)(b1::Baryon, b2::Baryon)
    b1.atoms == b2.atoms
end

function hash(b::Baryon, h::UInt)
    hash(b.atoms, h)
end

"""
    ismolecule(b::Baryon)

`true` when `b` is composed of more than one atom
"""
ismolecule(b::Baryon) = length(b.atoms) > 1

"""
    mass(b::Baryon)

Returns the mass [g] of `b`.
"""
function mass(b::Baryon)
    sum(atomic_masses[a] for a in b.atoms)
end

const roman_numerals = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"]
"""
Represents an atom or molecule (a `Baryon`) with a particular number of electrons (regardless of 
their configuration).
"""
struct Species
    baryon::Baryon
    charge::Int
end

"""
    Species(code::AbstractString)

Parse the "species code" as it is often specified in linelists and return a the "astronomy" 
notation. 01.00 → "H_I", 02.01 → "He_II", 02.1000 → "He_II", 0608 → "CO_I", etc.  
"""
function Species(code::AbstractString)
    code = strip(code, ['0', ' '])
    toks = split(code, [' ', '.', '_'])
    if length(toks) > 2
        throw(ArgumentError(code * " isn't a valid species code"))
    end
    baryon = Baryon(toks[1])
    charge = if length(toks) == 1 || length(toks[2]) == 0
        0 #no charge specified -> assume neutral
    else
        charge = findfirst(toks[2] .== roman_numerals)
        charge = (charge isa Int ? charge : parse(Int, toks[2]))
        #if this is a MOOG-style numeric code, the charge is correct, otherwise subtract 1
        if tryparse(Float64, code) == nothing 
            charge -= 1
        end
        charge
    end
    Species(baryon, charge)
end

function (==)(s1::Species, s2::Species)
    s1.baryon == s2.baryon && s1.charge == s2.charge
end

function hash(s::Species, h::UInt)
    hash((s.baryon, s.charge), h)
end

ismolecule(s::Species) = ismolecule(s.baryon)
mass(s::Species) = mass(s.baryon)

#This type represents an individual line.
struct Line{F} 
    wl::F                     #cm
    log_gf::F                 #unitless
    species::Species           
    E_lower::F                #eV (also called the excitation potential)
    gamma_rad::F              #s^-1
    gamma_stark::F            #s^-1
    vdW::Union{F, Tuple{F,F}} #either log(Γ_vdW) per electron or (σ, α) from ABO theory

    """
        Line(wl, log_gf, species, E_lower, gamma_rad, gamma_stark, vdW)

    Construct a `Line` with a possibly packed vdW parameter (sigma.alpha) format.  If vdW < 0,
    interpret it as log10(Γ) per particle.  Otherwise, interpret it as packed ABO parameters.
    """
    function Line(wl::F, log_gf::F, species::Species, E_lower::F, gamma_rad::F, gamma_stark::F,
                  vdW::F) where F <: Real
        new{F}(wl, log_gf, species, E_lower, gamma_rad, gamma_stark, 
               if vdW > 0
                   floor(vdW) * bohr_radius_cgs * bohr_radius_cgs, vdW - floor(vdW)
               elseif vdW == 0
                   0.0
               else 
                   10^vdW
               end
              )
    end
end
"""
    Line(wl, log_gf, species, E_lower)

Construct a `Line` without explicit broadening parameters.  They will be set automatically.
"""
function Line(wl::F, log_gf::F, species::Species, E_lower::F) where F <: Real
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
10,000 K. Used for atomic lines with no vdW and stark broadening info in the linelist.
Returns log10(γ_stark), log10(γ_vdW)

In the calculation of n*², uses the approximation that
\\overbar{r^2} = 5/2 {n^*}^4 / Z^2
which neglects the dependence on the angular momentum quantum number, l, in the the form given by
Warner 1967.

For autoionizing lines (those for which E_upper > χ), returns 0.0 for γ_vdW.
"""
function approximate_gammas(wl, species, E_lower; ionization_energies=ionization_energies)
    if ismolecule(species)
        return 0.0,0.0
    end

    Z = species.charge + 1 #Z is ionization stage, not atomic number
    χ = ionization_energies[species.baryon.atoms[1]][Z]
    c = c_cgs
    h = hplanck_eV
    k = kboltz_cgs
    E_upper = E_lower + (h * c / wl)

    nstar4_upper = (Z^2 * Rydberg_eV / (χ - E_upper))^2
    #From Cowley 1971
    γstark = 0.77e-18 * nstar4_upper * wl^2

    Δrbar2 = (5/2) * Rydberg_eV^2 * Z^2 * (1/(χ - E_upper)^2 - 1/(χ - E_lower)^2)
    if χ < E_upper
        γvdW = 0.0
    else
        #(log) γ_vdW From R J Rutten's course notes. An equivalent form can be found in Gray 2005.
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
    read_linelist(fname; format="kurucz")

Parse the provided linelist. in "Kurucz" format.
Pass `format="kurucz"` for a [Kurucz linelist](http://kurucz.harvard.edu/linelists.html),
`format="vald"` for a Vald linelist, and `format="moog"` for a MOOG linelist.

Note that dissociation energies in a MOOG linelist will be ignored.
"""
function read_linelist(fname::String; format="vald") :: Vector{Line}
    format = lowercase(format)
    linelist = open(fname) do f
        if format == "kurucz"
            parse_kurucz_linelist(f)
        elseif format == "vald"
            parse_vald_linelist(f)
        elseif format == "moog"
            parse_moog_linelist(f)
        else
            throw(ArgumentError("$(format) is not a supported linelist format"))
        end
    end

    filter!(linelist) do line
        0 <= line.species.charge <= 2
    end

    #ensure linelist is sorted
    if !issorted(linelist, by=l->l.wl)
        sort!(linelist, by=l->l.wl)
    end

    linelist
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
            Species(line[19:24]),
            min(E_levels...),
            expOrZero(parse(Float64, line[81:86])),
            expOrZero(parse(Float64, line[87:92])),
            parse(Float64, line[93:98]))
    end
end

function parse_vald_linelist(f)
    lines = collect(eachline(f))

    #figure out how big the header is
    firstline = findfirst(lines) do line
        line[1] == '\''
    end
    header = split(lines[firstline-1])

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

    #air or vacuum wls?
    if contains(header[3], "air")
        wl_transform = air_to_vacuum
    elseif contains(header[3], "vac")
        wl_transform = identity
    else
        throw(ArgumentError(
            "Can't parse linelist.  I don't understand this wavelength column name: " * header[3]))
    end
    #Energy in cm^-1 or eV?
    E_col = header[shortformat ? 4 : 6]
    if contains(E_col, "eV")
        E_transform = identity
    elseif contains(E_col, "cm")
        E_transform(x) = x * c_cgs * hplanck_eV
    else
        throw(ArgumentError(
            "Can't parse linelist.  I don't understand this energy column name: " * E_col))
    end

    map(lines) do line
        toks = split(line, ',')
        if shortformat
            #extract all
            if firstline == 3
                new_line_imputing_zeros(
                     wl_transform(parse(Float64, toks[2])*1e-8),
                     parse(Float64, toks[4]),
                     Species(strip(toks[1], ['\''])),
                     E_transform(parse(Float64, toks[3])),
                     expOrZero(parse(Float64, toks[5])),
                     expOrZero(parse(Float64, toks[6])),
                     parse(Float64, toks[7]))
            #extract stellar
            elseif firstline == 4
                new_line_imputing_zeros(
                     wl_transform(parse(Float64, toks[2])*1e-8),
                     parse(Float64, toks[5]),
                     Species(strip(toks[1], ['\''])),
                     E_transform(parse(Float64, toks[3])),
                     expOrZero(parse(Float64, toks[6])),
                     expOrZero(parse(Float64, toks[7])),
                     parse(Float64, toks[8]))
            else
                throw(ArgumentError("Can't determine if this is an \"extract all\" or \"extract " *
                                    "stellar\" format linelist"))
            end
        else
            new_line_imputing_zeros(
                parse(Float64, toks[2])*1e-8,
                parse(Float64, toks[3]),
                Species(strip(toks[1], ['\''])),
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
        Line(parse(Float64, toks[1]) * 1e-8, #convert Å to cm
             parse(Float64, toks[4]),
             Species(toks[2]),
             parse(Float64, toks[3]))
    end
    #TODO issue warning, don't autoconvert
    linelist
end
