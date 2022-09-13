using StaticArrays
using CSV

"""
Represents an atom or molecule, irespective of its charge.
"""
struct Formula
    #supports up to triatomic molecules, can be trivially extended.
    #Unlike tuples SVectors support sorting
    atoms::SVector{3, UInt8}

    function Formula(Z::Integer) 
        @assert 1 <= Z <= Natoms
        new([0x00, 0x00, Z])
    end
    function Formula(Zs::AbstractVector{<:Integer})
        l = length(Zs)
        if l == 0
            throw(ArgumentError("Can't construct an empty Formula"))
        elseif l == 1
            Zs = [0 ; 0 ; Zs]
        elseif l == 2 
            Zs = [0 ; Zs]
        elseif l > 3
            throw(ArgumentError("Can't construct Formula with more than three atoms"))
        end
        @assert(issorted(Zs))
        new(Zs)
    end

    """
        Formula(code::String)

        Construct a Formula from an encoded string form.  This can be a MOOG-style numeric code, i.e.
        "0801" for OH, or an atomic or molecular symbol, i.e. "FeH", "Li", or "C2".
    """
    function Formula(code::String) 
        if code in atomic_symbols #quick-parse single elements
            return new([0, 0, atomic_numbers[code]]) 
        end

        #handle numeric codes, e.g. 0801 -> OH
        if all(isdigit(c) for c in code)
            if length(code) <= 2
                return Formula(parse(Int,code))
            elseif length(code) <= 4
                el1 = parse(Int, code[1:end-2])   #first digit
                el2 = parse(Int, code[end-1:end]) #second digit
                return new([0x00, min(el1, el2), max(el1, el2)])
            else
                throw(ArgumentError("numeric codes for molecules with more than 4 chars like " * 
                                    "$(code) are not supported"))
            end
        end

        #otherwise, code should be "OH", "FeH", "Li", "C2", etc.
        inds::Vector{Int} = findall(code) do c
            isdigit(c) || isuppercase(c)
        end
        push!(inds, length(code)+1)
        subcode::Vector{String} = map(1:(length(inds)-1)) do j
            code[inds[j]:inds[j+1]-1]
        end

        atoms = UInt8[]
        for s in subcode
            num = tryparse(Int, s)
            if num isa Int
                previous = atoms[end]
                for _ in 1:(num-1)
                    push!(atoms, previous)
                end
            else
                push!(atoms, atomic_numbers[s])
            end
        end

        sort!(atoms)
        Formula(atoms)
    end
end

function get_atoms(f::Formula) 
    if f.atoms[2] == 0
        view(f.atoms,3:3)
    elseif f.atoms[1] == 0
        view(f.atoms,2:3)
    else
        view(f.atoms, 1:3)
    end
end

# it's important that this produces something parsable by the constructor
function Base.show(io::IO, f::Formula)
    print(io, *([atomic_symbols[i] for i in f.atoms if i != 0]...))
end

"""
    ismolecule(f::Formula)

`true` when `f` is composed of more than one atom
"""
ismolecule(f::Formula) = f.atoms[2] != 0

"""
    get_mass(f::Formula)

Returns the mass [g] of `f`.
"""
function get_mass(f::Formula)
    sum(atomic_masses[a] for a in get_atoms(f))
end

const roman_numerals = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"]
"""
Represents an atom or molecule (a `Formula`) with a particular number of electrons (regardless of 
their configuration).
"""
struct Species
    formula::Formula
    charge::Int
end

"""
    Species(code::AbstractString)

Parse the "species code" as it is often specified in linelists and return a the "astronomy" 
notation. 01.00 → "H I", 02.01 → "He II", 02.1000 → "He II", 0608 → "CO I", etc.  

TODO


To parse at compile time, use the `species` string macro, i.e. `species"H I"`.
"""
function Species(code::AbstractString)
    code = strip(code, ['0', ' '])
    toks = split(code, [' ', '.', '_'])
    filter!(toks) do tok
        tok != "" #TODO add tests
    end
    if length(toks) > 2
        throw(ArgumentError(code * " isn't a valid species code"))
    end
    #convert toks[1] from Substring to String.  Better for type stability in Formula
    formula = Formula(String(toks[1])) 
    charge = if length(toks) == 1 || length(toks[2]) == 0
        0 #no charge specified -> assume neutral
    else
        charge = findfirst(toks[2] .== roman_numerals)
        charge = (charge isa Int ? charge : parse(Int, toks[2]))
        #if this is a MOOG-style numeric code, the charge is correct, otherwise subtract 1
        if tryparse(Float64, code) === nothing 
            charge -= 1
        end
        charge
    end
    Species(formula, charge)
end

#used to contruct Species at compile time and avoid parsing in hot loops
macro species_str(s)
    Species(s)
end

function Base.show(io::IO, s::Species)
    show(io, s.formula)
    print(io, " ", get_roman_numeral(s))
end

ismolecule(s::Species) = ismolecule(s.formula)
get_mass(s::Species) = get_mass(s.formula)
get_atoms(s::Species) = get_atoms(s.formula)
get_roman_numeral(s::Species) = get(roman_numerals,s.charge+1, string(s.charge+1))

#This type represents an individual line.
struct Line{F} 
    wl::F                     #cm
    log_gf::F                 #unitless
    species::Species           
    E_lower::F                #eV (also called the excitation potential)
    gamma_rad::F              #s^-1
    gamma_stark::F            #s^-1
    vdW::Union{F, Tuple{F,F}} #either γ_vdW [s^-1] per electron or (σ, α) from ABO theory

    @doc """
        Line(wl::F, log_gf::F, species::Species, E_lower::F, 
             gamma_rad::Union{F, Missing}=missing, gamma_stark::Union{F, Missing}=missing, 
             vdw::Union{F, Tuple{F, F}, Missing}, missing) where F <: Real

    Arguments:
     - `wl`: wavelength, in cm
     - `log_gf`: (log base 10) oscillator strength (unitless)
     - `species`: the `Species` associated with the line
     - `E_lower`: The energy (excitiation potential) of the lower energy level (eV)

    Optional Arguments (these override default recipes):
     - `gamma_rad`: Fundemental width
     - `gamma_stark`: Stark broadening width at 10,000 K (s⁻¹)
     - `vdW`: Either the van der Waals broadening width at 10,000 K (s⁻¹) or a `Tuple`, (σ, α) from
       ABO theory.

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

# it's important that this produces something parsable by the constructor
function Base.show(io::IO, line::Line)
    show(io, line.species)
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
    approximate_gammas(wl, species, E_lower; ionization_energies=Korg.ionization_energies)

A simplified form of the Unsoeld (1955) approximation for van der Waals broadening and the 
[Cowley 1971](https://ui.adsabs.harvard.edu/abs/1971Obs....91..139C/abstract) approximation for 
Stark broadening, evaluated at 10,000 K. 
Used for atomic lines with no vdW and stark broadening info in the linelist.

Returns `(γ_stark`, `log10(γ_vdW))` in Hz, where these are the per-perturber quantities.
For autoionizing lines (those for which E_upper > χ), Returns 0.0 for γ_vdW.

In the calculation of `n*²`, uses the approximation that
``\\overbar{r^2} = 5/2 {n^*}^4 / Z^2``
which neglects the dependence on the angular momentum quantum number, l, in the the form given by
[Warner 1967](https://ui.adsabs.harvard.edu/abs/1967MNRAS.136..381W/abstract) (the earliest english 
work reporting the Unsoeld result).
"""
function approximate_gammas(wl, species, E_lower; ionization_energies=ionization_energies)
    Z = species.charge + 1 #Z is ionization stage, not atomic number
    if ismolecule(species) || Z > 3
        return 0.0,0.0
    end
    χ = ionization_energies[get_atoms(species.formula)[1]][Z]
    c = c_cgs
    h = hplanck_eV
    k = kboltz_cgs
    E_upper = E_lower + (h * c / wl)

    #It's not obvious to me which Rydberg constant to use here, and below in Δrbar2.  The sources
    #are not entirely clear. It doesn't make a big difference.
    nstar4_upper = (Z^2 * RydbergH_eV / (χ - E_upper))^2
    #I'm not actually able to reproduce Crowley 1971 equation 7 (his simplified form) from equation 
    #5, but these match the values in the Turbospectrum source, so they are probably correct.
    #The constants here were calculated assuming that "v" is the mean (not modal) electron speed
    if Z == 1
        γstark = 2.25910152e-7 * nstar4_upper #Cowley (1971) equation 5 evaluated at T=10,000 K
    else
        #Cowley (1971) equation 6 @ T=10,000 K (n.b. the constant is 12/5 * that above)
        γstark = 5.42184365e-7 * nstar4_upper / (Z + 1)^2 
    end

    Δrbar2 = (5/2) * Rydberg_eV^2 * Z^2 * (1/(χ - E_upper)^2 - 1/(χ - E_lower)^2)
    if χ < E_upper
        γvdW = 0.0
    else
        # (log) γ_vdW From R J Rutten's course notes. 
        # Equations 11.29 and 11.30 from Gray 2005 are equivalent 
        γvdW = 6.33 + 0.4log10(Δrbar2) + 0.3log10(10_000) + log10(k)
    end

    γstark, γvdW
end

"""
    read_linelist(filename; format="vald", isotopic_abundances=Korg.isotopic_abundances)

Parse a linelist file, returning a vector of [`Line`](@ref)s.

Pass `format="kurucz"` for a [Kurucz linelist](http://kurucz.harvard.edu/linelists.html) 
(`format=kurucz_vac` if it uses vacuum wavelengths; Be warned that Korg will not assume that 
wavelengths are vacuum below 2000 Å),`format="vald"` for a 
[VALD](http://vald.astro.uu.se/~vald/php/vald.php) linelist, and `format="moog"` for a MOOG linelist
(doesn't support broadening parameters or dissociation energies).  

VALD linelists (the default and preferred format) can be either "short" or "long" format, 
"extract all" or "extract stellar".  Air wavelengths will automatically be converted into vacuum
wavelengths, and energy levels will be automatically converted from cm``^{-1}`` to eV.

When they are not pre-scaled by isotopic abundace (which VALD does by default), Korg will 
automatically adjust the `log_gf` of each line according to `isotopic_abundances`, which defaults 
to the values from [NIST](https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses).
To use custom isotopic abundances, just pass `isotopic_abundances` as a dictionary mapping
`(atomic number, atomic weight)` pairs to abundances between 0 and 1.
"""
function read_linelist(fname::String; format="vald", isotopic_abundances=isotopic_abundances)
    format = lowercase(format)
    linelist = open(fname) do f
        if format == "kurucz"
            parse_kurucz_linelist(f; vacuum=false)
        elseif format == "kurucz_vac"
            parse_kurucz_linelist(f; vacuum=true)
        elseif format == "vald"
            parse_vald_linelist(f, isotopic_abundances)
        elseif format == "moog"
            parse_moog_linelist(f)
        else
            throw(ArgumentError("$(format) is not a supported linelist format"))
        end
    end

    filter!(linelist) do line #filter triply+ ionized and hydrogen lines
        (0 <= line.species.charge <= 2) && (line.species != species"H_I")
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

function parse_kurucz_linelist(f; vacuum=false)
    lines = Line[]
    for row in eachline(f)
        row == "" && continue #skip empty lines

        #some linelists have a missing column in the wavelenth region
        if length(row) == 159 
            row = " " * row
        end
        
        #kurucz provides wavenumbers for "level 1" and "level 2", which is which is 
        #determined by parity
        E_levels = map((row[25:36], row[53:64])) do s
            #abs because Kurucz multiplies predicted values by -1
            abs(parse(Float64,s)) * c_cgs * hplanck_eV
        end

        wl_transform = vacuum ? identity : air_to_vacuum

        push!(lines, Line(wl_transform(parse(Float64, row[1:11])*1e-7), #convert from nm to cm
                     parse(Float64, row[12:18]),
                     Species(row[19:24]),
                     min(E_levels...),
                     expOrMissing(parse(Float64, row[81:86])),
                     expOrMissing(parse(Float64, row[87:92])),
                     idOrMissing(parse(Float64, row[93:98]))))
    end
    lines
end

function parse_vald_linelist(f, isotopic_abundances)
    lines = filter!(collect(eachline(f))) do line
        length(line) > 0 && line[1] != '#' #remove comments and empty lines
    end

    lines = replace.(lines, "'"=>"\"") #replace single quotes with double

    # is this an "extract all" or an "extract stellar" linelist?
    extractall = !occursin(r"^\s+\d", lines[1])
    firstline = extractall ? 3 : 4
    header = lines[firstline - 1]

    scale_isotopes = any(startswith.(lines, "* oscillator strengths were NOT scaled "))
    if !scale_isotopes && !any(startswith.(lines,"* oscillator strengths were scaled "))
        throw(ArgumentError("Can't parse linelist.  Can't detect whether log(gf)s are scaled by "*
                            "isotopic abundance."))
    end

    #we take the linelist to be long-format when the second line after the header starts with a 
    #space or a quote followed a space.  In some linelists the quotes are there, but in others 
    #they are not.
    shortformat = !(occursin(r"^\"? ", lines[firstline + 1])) 
    body = lines[firstline : (shortformat ? 1 : 4) : end]
    body = body[1 : findfirst(l->l[1]!='\"' || !isuppercase(l[2]), body)-1]

    CSVheader = if shortformat && extractall
        ["species", "wl", "E_low", "loggf", "gamma_rad", "gamma_stark", "gamma_vdW", "lande", 
         "reference"]
    elseif shortformat #extract stellar
        ["species", "wl", "E_low", "Vmic", "loggf", "gamma_rad", "gamma_stark", "gamma_vdW", 
         "lande", "depth", "reference"]
    else #long format (extract all or extract stellar)
        ["species", "wl", "loggf", "E_low", "J_lo", "E_up", "J_up", "lower_lande", "upper_lande",
         "mean_lande", "gamma_rad", "gamma_stark", "gamma_vdW"]
    end
    body = CSV.File(reduce(vcat, codeunits.(body.*"\n")), header=CSVheader, delim=',', 
                    silencewarnings=true)

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

    loggf = body.loggf .+ if scale_isotopes
        refs = if !shortformat #the references are on different lines
            lines[firstline+3 .+ ((0:length(body)-1) .* 4)]
        else #references are in the last column
            body.reference
        end

        map(refs) do ref
            #find things that look like (16)O or (64)Ni in reference string
            regexp = r"\((?<isotope>\d\d?\d?)\)(?<elem>\p{Lu}\p{Ll}?)"
            #add up the adjustments to log(gf) from isotopic abundances (zero if no info is present)
            log_probs = map(findall(regexp, ref)) do r
                m = match(regexp, ref[r])
                log10(isotopic_abundances[atomic_numbers[m["elem"]], parse(Int64, m["isotope"])])
            end
            sum([0 ; log_probs])
        end
    else
        0
    end

    Line.(wl * 1e-8, loggf, Species.(body.species), E_low, expOrMissing.(body.gamma_rad), 
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
