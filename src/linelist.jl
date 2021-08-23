import Base: (==), hash
using CSV

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
        if code in atomic_symbols
            return new([code]) #quick-parse single elements
        end

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

    function Baryon(Z::Int...)
        new([atomic_symbols[z] for z in sort(collect(Z))])
    end

    function Baryon(symbols::Vector{AbstractString})
        @assert issorted(symbols, by=s->atomic_numbers[s])
        new(String.(symbols))
    end
end

#pretty-print lines in REPL and jupyter notebooks
function Base.show(io::IO, m::MIME"text/plain", b::Baryon)
    print(io, *(b.atoms...))
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

#useful constants
const H_I = Species("H_I")
const H_II = Species("H_II")
const He_I = Species("He_I")
const He_II = Species("He_II")
const He_III = Species("He_III")

#pretty-print lines in REPL and jupyter notebooks
function Base.show(io::IO, m::MIME"text/plain", s::Species)
    show(io, m, s.baryon)
    print(io, " ", roman_numerals[s.charge+1])
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
    vdW::Union{F, Tuple{F,F}} #either Γ_vdW [s^-1] per electron or (σ, α) from ABO theory

    """
        Line(wl::F, log_gf::F, species::Species, E_lower::F, 
             gamma_rad::Union{F, Missing}=missing, gamma_stark::Union{F, Missing}=missing, 
             vdw::Union{F, Tuple{F, F}, Missing}, missing) where F <: Real

    Construct a `Line`.  If any of `gamma_rad`, `gamma_stark`, or `vdW` are `missing`, guess them.
    `vdW` may be log(Γ_vdW) (assumed if negative), Γ_vdW (assumed if 0 < `vdW` < 1), or packed ABO 
    parameters (assumed if `vdW` > 1).  It may also be passed as a Tuple, `(σ, α)`.
    """
    function Line(wl::F, log_gf::F, species::Species, E_lower::F, 
                  gamma_rad::Union{F, Missing}=missing, gamma_stark::Union{F, Missing}=missing, 
                  vdW::Union{F, Tuple{F, F}, Missing}=missing) where F <: Real
        if ismissing(gamma_stark) || ismissing(vdW)
            gamma_stark_approx, vdW_approx = approximate_gammas(wl, species, E_lower)
            if ismissing(gamma_stark)
                gamma_stark = gamma_stark_approx
            end
            if ismissing(vdW)
                vdW = vdW_approx
            end
        end
        if ismissing(gamma_rad)
            gamma_rad = approximate_radiative_gamma(wl, log_gf)
        end
        
        if vdW isa F
            if vdW < 0 #if vdW is negative, assume it's log(Γ_vdW) 
                vdW = 10^vdW
            elseif vdW > 1 #if it's > 1 assume it's packed ABO params
                vdW = (floor(vdW) * bohr_radius_cgs * bohr_radius_cgs, vdW - floor(vdW))
            end
        end 

        new{F}(wl, log_gf, species, E_lower, gamma_rad, gamma_stark, vdW)
    end
end

#pretty-print lines in REPL and jupyter notebooks
function Base.show(io::IO, m::MIME"text/plain", line::Line)
    show(io, m, line.species)
    print(io, " ", round(line.wl*1e8, digits=6), " Å")
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
    Z = species.charge + 1 #Z is ionization stage, not atomic number
    if ismolecule(species) || Z > 3
        return 0.0,0.0
    end
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

    filter!(linelist) do line #filter triply+ ionized and hydrogen lines
        (0 <= line.species.charge <= 2) && (line.species.baryon.atoms != ["H"])
    end

    #ensure linelist is sorted
    if !issorted(linelist, by=l->l.wl)
        sort!(linelist, by=l->l.wl)
    end

    linelist
end

#used to handle missing gammas in vald and kurucz lineslist parsers
expOrMissing(x) = x == 0.0 ? missing : 10.0^x
idOrMissing(x) = x == 0.0 ? missing : x

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
             Species(line[19:24]),
             min(E_levels...),
             expOrMissing(parse(Float64, line[81:86])),
             expOrMissing(parse(Float64, line[87:92])),
             idOrMissing(parse(Float64, line[93:98])))
    end
end

function parse_vald_linelist(f)
    lines = collect(eachline(f))

    # is this an "extract all" or an "extract stellar" linelist?
    extractall = !occursin(r"^\s+\d", lines[1])
    firstline = extractall ? 3 : 4
    header = lines[firstline - 1]

    isotope_scaled = if "* oscillator strengths were scaled by the solar isotopic ratios." in lines
        true
    elseif "* oscillator strengths were NOT scaled by the solar isotopic ratios." in lines
        throw("Isotopic scaling not yet implemented.")
    else
        throw(ArgumentError("Can't parse linelist.  Can't detect whether log(gf)s are scaled by "*
                            ":isotopic abundance."))
    end

    #vald short or long format?
    shortformat = isuppercase(lines[firstline][2]) && isuppercase(lines[firstline+1][2])

    body = lines[firstline : (shortformat ? 1 : 4) : end]
    body = body[1 : findfirst(l->l[1]!='\'' || !isuppercase(l[2]), body)-1]
    
    CSVheader = if shortformat && extractall
        ["species","wl","E_low","loggf","gamma_rad","gamma_stark","gamma_vdW"]
    elseif shortformat #extract stellar
        ["species","wl", "E_low","Vmic","loggf","gamma_rad","gamma_stark","gamma_vdW"]
    elseif extractall #long format
        ["species","wl","loggf","E_low","E_up","lande","gamma_rad","gamma_stark","gamma_vdW"]
    else #long format extract stellar
        ["species","wl","loggf","E_low","J_lo","E_up","J_up","lower_lande","upper_lande",
         "mean_lande","gamma_rad","gamma_stark","gamma_vdW"]
    end
    body = CSV.File(reduce(vcat, codeunits.(body.*"\n")), header=CSVheader, silencewarnings=true)

    species = (s->s[2:end-1]).(body.species) #strip quotes

    E_low = if contains(header, "cm") #convert E_low to eV if necessary
        body.E_low * c_cgs * hplanck_eV
    elseif contains(header, "eV")
        body.E_low
    else
        throw(ArgumentError( "Can't parse linelist.  Can't determine energy units: " * E_col))
    end

    wl = if contains(header, "air") #convert wls to vacuum if necessary
        air_to_vacuum.(body.wl)
    elseif contains(header, "vac")
        body.wl
    else
        throw(ArgumentError( "Can't parse linelist.  Can't determine vac/air wls: " * header))
    end

    Line.(wl * 1e-8, body.loggf, Species.(species), E_low, expOrMissing.(body.gamma_rad), 
          expOrMissing.(body.gamma_stark), idOrMissing.(body.gamma_vdW))
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
    linelist
end
