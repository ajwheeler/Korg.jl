using StaticArrays

const MAX_ATOMS_PER_MOLECULE = 6

"""
Represents an atom or molecule, irrespective of its charge.
"""
struct Formula
    # Support molecules with up to MAX_ATOMS_PER_MOLECULE atoms. 
    # Unlike tuples, SVectors support sorting, which is why we use them here.
    atoms::SVector{6,UInt8}

    function Formula(Z::Integer)
        @assert 1 <= Z <= MAX_ATOMIC_NUMBER
        new([zeros(UInt8, MAX_ATOMS_PER_MOLECULE - 1); Z])
    end
    function Formula(Zs::AbstractVector{<:Integer})
        len = length(Zs)
        if len == 0
            throw(ArgumentError("Can't construct an empty Formula"))
        elseif len < MAX_ATOMS_PER_MOLECULE
            Zs = [zeros(Int, MAX_ATOMS_PER_MOLECULE - len); Zs]
        elseif len > MAX_ATOMS_PER_MOLECULE
            throw(ArgumentError("Can't construct Formula with atoms $(Int.(Zs)). Up to $(MAX_ATOMS_PER_MOLECULE) atoms are supported."))
        end
        @assert(issorted(Zs))
        new(Zs)
    end

    """
        Formula(code::String)

    Construct a Formula from an encoded string form.  This can be a MOOG-style numeric code, 
    i.e. "0801" for OH, or an atomic or molecular symbol, i.e. "FeH", "Li", or "C2".
    """
    function Formula(code::AbstractString)
        if code in atomic_symbols #quick-parse single elements
            return new([zeros(UInt8, MAX_ATOMS_PER_MOLECULE - 1); atomic_numbers[code]])
        end

        #handle numeric codes, e.g. 0801 -> OH
        if all(isdigit(c) for c in code)
            if length(code) <= 2
                Formula(parse(Int, code))
            elseif length(code) <= 4
                el1 = parse(Int, code[1:end-2])   #first digit
                el2 = parse(Int, code[end-1:end]) #second digit
                new([zeros(UInt8, MAX_ATOMS_PER_MOLECULE - 2); min(el1, el2); max(el1, el2)])
            elseif length(code) / 2 <= MAX_ATOMS_PER_MOLECULE
                # if there are an odd number of digits, pad the front with a zero
                if length(code) % 2 == 1
                    code = "0" * code
                end
                els = map(1:2:(length(code)-1)) do i
                    parse(Int, code[i:i+1])
                end
                new([zeros(UInt8, MAX_ATOMS_PER_MOLECULE - length(els)); sort(els)])
            else
                throw(ArgumentError("Korg only supports atoms with up to $MAX_ATOMS_PER_MOLECULE nuclei. (Trying to parse $code)"))
            end
        else
            #otherwise, code should be "OH", "FeH", "Li", "C2", etc.
            inds::Vector{Int} = findall(code) do c
                isdigit(c) || isuppercase(c)
            end
            push!(inds, length(code) + 1)
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
end

"""
    get_atoms(x)

Returns an array view containing the atomic number of each atom that makes up the formula or species
x.  E.g. `get_atoms(Korg.species"H2O")` yields [1, 1, 8].
"""
function get_atoms(f::Formula)
    i = findlast(f.atoms .== 0)
    if isnothing(i)
        view(f.atoms, 1:MAX_ATOMS_PER_MOLECULE)
    else
        view(f.atoms, i+1:MAX_ATOMS_PER_MOLECULE)
    end
end

"""
    get_atom(x)

Returns the atomic number of an atomic Korg.Species or Korg.Formula.
"""
function get_atom(f::Formula)
    if ismolecule(f)
        throw(ArgumentError("Can't get the atomic number of a molecule.  Use `get_atoms` instead."))
    end
    get_atoms(f)[1]
end

"""
    n_atoms(x)

The number of atoms in the Korg.Species or Korg.Formula x.
"""
function n_atoms(f::Formula)
    i = findlast(f.atoms .== 0)
    if isnothing(i)
        MAX_ATOMS_PER_MOLECULE
    else
        MAX_ATOMS_PER_MOLECULE - i
    end
end

# it's important that this produces something parsable by the constructor
Base.show(io::IO, f::Formula) = print(io, *([atomic_symbols[i] for i in f.atoms if i != 0]...))

# make it broadcast like a scalar
Base.broadcastable(f::Formula) = Ref(f)

"""
    ismolecule(f::Formula)

`true` when `f` is composed of more than one atom
"""
ismolecule(f::Formula) = f.atoms[MAX_ATOMS_PER_MOLECULE-1] != 0

"""
    get_mass(f::Formula)

Returns the mass [g] of `f`.
"""
get_mass(f::Formula) = sum(atomic_masses[a] for a in get_atoms(f))

"""
Represents an atom or molecule (a `Formula`) with a particular number of electrons (regardless of
their configuration).
"""
struct Species
    formula::Formula
    charge::Int16

    function Species(f::Formula, charge::Integer)
        if charge < -1
            throw(ArgumentError("Can't construct a species with charge < -1: $(f) with charge $charge"))
        end
        new(f, Int16(charge))
    end
end

"""
    Species(code)

Parse the "species code" in many of the forms in which it is often specified and return an object
representing the species. `code` can be either a string or a float.

# Examples

  - "H I" -> H I
  - "H 1" -> H I
  - "H     1" -> H I
  - "H_1" -> H I
  - "H.I" -> H I
  - "H 2" -> H II
  - "H2" -> H₂
  - "H" -> H I
  - "01.00" → H I
  - "02.01" → He II
  - "02.1000" → He II
  - "0608" → CO I

!!! note

    To parse at compile time, use the `species` string macro, i.e. `species"H I"`.  This is
    important in hot inner loops.

!!! warning

    MOOG codes which include isotopic information will not be parsed correctly by this function,
    though [`read_linelist`](@ref) handles them correctly. # leading 0s are safe to remove
"""
function Species(code::AbstractString)
    code = strip(code, ['0', ' ']) # leading 0s are safe to remove

    # if the species ends in "+" or "-", convert it to a numerical charge. Remember, the ionization 
    # number is the charge+1, so for us "H 0" is H⁻ and "H 2" in H⁺.
    if code[end] == '+'
        code = code[1:end-1] * " 2"
    elseif code[end] == '-'
        code = code[1:end-1] * " 0"
    end

    # these are the valid separators between the atomic number part of a species code and the 
    # charge-containing part.  For example, "01.01" parses the same as "01 01", but "01,01" fails.
    # Or, "C 2" parses the same at "C.2".  "-"s can't be separators as they can be minus signs.
    toks = split(code, [' ', '.', '_'])
    # this allows for leading, trailing, and repeat separators.  "01.01" parses the same as 
    # ".01..01.".
    filter!(!=(""), toks)

    if length(toks) > 2
        throw(ArgumentError(code * " isn't a valid species code"))
    end
    #convert toks[1] from Substring to String.  Better for type stability in Formula
    formula = Formula(String(toks[1]))
    charge = if length(toks) == 1 || length(toks[2]) == 0
        0 #no charge specified -> assume neutral
    else
        # first check if the "charge tag" is a roman numeral.  If it's not, parse it as an Int.
        charge = findfirst(toks[2] .== roman_numerals)
        charge = (charge isa Int ? charge : parse(Int, toks[2]))
        # if this is a Kurucz-style numeric code, the charge is correct, otherwise subtract 1
        if tryparse(Float64, code) === nothing
            charge -= 1
        end
        charge
    end
    Species(formula, charge)
end
Species(code::AbstractFloat) = Species(string(code))

#used to contruct Species at compile time and avoid parsing in hot loops
macro species_str(s)
    Species(s)
end

const roman_numerals = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"]
function Base.show(io::IO, s::Species)
    show(io, s.formula)
    if ismolecule(s) && s.charge == 1
        print(io, "+")
    elseif ismolecule(s) && s.charge == 0
        # no charge tag for neutral molecules
    elseif 0 <= s.charge <= length(roman_numerals) - 1
        print(io, " ", roman_numerals[s.charge+1])
    elseif s.charge == -1
        print(io, "-")
    else
        print(io, " ", s.charge)
    end
end

# make it broadcast like a scalar
Base.broadcastable(s::Species) = Ref(s)

ismolecule(s::Species) = ismolecule(s.formula)
get_mass(s::Species) = get_mass(s.formula)
get_atoms(s::Species) = get_atoms(s.formula)
get_atom(s::Species) = get_atom(s.formula)
n_atoms(s::Species) = n_atoms(s.formula)

"""
    all_atomic_species()

Returns an iterator that runs over all atomic species supported by Korg.
"""
all_atomic_species() = (Species(Formula(Z), charge)
                        for Z in 1:MAX_ATOMIC_NUMBER, charge in 0:2 if charge <= Z)
