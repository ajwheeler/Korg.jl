using CSV, HDF5, LazyArtifacts
using Pkg.Artifacts: @artifact_str

#This type represents an individual line.
struct Line{F1,F2,F3,F4,F5,F6}
    wl::F1                       #cm
    log_gf::F2                   #unitless
    species::Species
    E_lower::F3                  #eV (also called the excitation potential)
    gamma_rad::F4                #s^-1
    gamma_stark::F5              #s^-1
    vdW::Union{F6,Tuple{F6,F6}} #either γ_vdW [s^-1] per electron or (σ, α) from ABO theory

    @doc """
        Line(wl::F, log_gf::F, species::Species, E_lower::F, 
             gamma_rad::Union{F, Missing}=missing, gamma_stark::Union{F, Missing}=missing, 
             vdw::Union{F, Tuple{F, F}, Missing}, missing) where F <: Real

    Arguments:
     - `wl`: wavelength (Assumed to be in cm if < 1, otherwise in Å)
     - `log_gf`: (log base 10) oscillator strength (unitless)
     - `species`: the `Species` associated with the line
     - `E_lower`: The energy (excitiation potential) of the lower energy level (eV)

    Optional Arguments (these override default recipes):
     - `gamma_rad`: Fundemental width
     - `gamma_stark`: per-perturber Stark broadening width at 10,000 K (s⁻¹).
     - `vdW`: If this is present, it may may be 
         - `log10(γ_vdW)`, assumed if negative
         - 0, corresponding to no vdW broadening
         - A fudge factor for the Unsoeld approximation, assumed if between 0 and 20
         - The [ABO parameters](https://www.astro.uu.se/~barklem/howto.html) as packed float 
           (assumed if >= 20) or a `Tuple`, `(σ, α)`.

        This behavior is intended to mirror that of Turbospectrum as closely as possible.


    See [`approximate_gammas`](@ref) for more information on the default recipes for `gamma_stark` 
    and `vdW`.

    Note the the "gamma" values here are FWHM, not HWHM, of the Lorenztian component of the line 
    profile, and are in units of s⁻¹.


        Line(line::line; kwargs...)

    Construct a new `Line` by copying the values from an existing `Line`.  Any of the values can be 
    modified with keyword arguments, e.g. `Line(line, log_gf=0.0)`.
    """
    function Line(wl::F1, log_gf::F2, species::Species, E_lower::F3,
                  gamma_rad::Union{F4,Missing}=missing, gamma_stark::Union{F5,Missing}=missing,
                  vdW::Union{F6,Tuple{F6,F6},Missing}=missing) where {F1<:Real,F2<:Real,F3<:Real,
                                                                      F4<:Real,F5<:Real,F6<:Real}
        if wl >= 1
            wl *= 1e-8 #convert Å to cm
        end
        if ismissing(gamma_stark) || isnan(gamma_stark) || ismissing(vdW) ||
           (!(vdW isa Tuple) && isnan(vdW))
            gamma_stark_approx, vdW_approx = approximate_gammas(wl, species, E_lower)
            if ismissing(gamma_stark) || isnan(gamma_stark)
                gamma_stark = gamma_stark_approx
            end
            if ismissing(vdW) || isnan(vdW)
                vdW = vdW_approx
            end
        end
        if ismissing(gamma_rad) || isnan(gamma_rad)
            gamma_rad = approximate_radiative_gamma(wl, log_gf)
        end

        # if vdW is a tuple, assume it's (σ, α) from ABO theory
        # if it's a float, there are four possibilities
        if !(vdW isa Tuple) #F6 will not be defined if vdW is missing
            vdW = if vdW < 0
                10^vdW  # if vdW is negative, assume it's log(γ_vdW) 
            elseif vdW == 0
                0.0  # if it's exactly 0, leave it as 0 (no vdW broadening)
            elseif vdW < 1e-2
                # if it's between 0 and 1e-2, assume it's γ_vdW
                vdW
            elseif vdW < 20
                # if it's between 0 and 20, assume it's a fudge factor for the Unsoeld approximation
                vdW * 10^(approximate_gammas(wl, species, E_lower)[2])
            else #if it's >= 20 assume it's packed ABO params
                (floor(vdW) * bohr_radius_cgs * bohr_radius_cgs, vdW - floor(vdW))
            end
        end

        new{F1,F2,F3,typeof(gamma_rad),typeof(gamma_stark),eltype(vdW)}(wl, log_gf, species,
                                                                        E_lower, gamma_rad,
                                                                        gamma_stark, vdW)
    end
end
# constructor to allow for copying a line and modifying some values (see docstring)
function Line(line::Line; wl=line.wl, log_gf=line.log_gf, species=line.species,
              E_lower=line.E_lower, gamma_rad=line.gamma_rad, gamma_stark=line.gamma_stark,
              vdW=line.vdW)
    Line(wl, log_gf, species, E_lower, gamma_rad, gamma_stark, vdW)
end

function Base.show(io::IO, ::MIME"text/plain", line::Line)
    show(io, line.species)
    print(io, " ", round(line.wl * 1e8; digits=6), " Å (log gf = ", round(line.log_gf; digits=2),
          ")")
end

# make it broadcast like a scalar
Base.broadcastable(l::Line) = Ref(l)

"""
    approximate_radiative_gamma(wl, log_gf)

Approximate radiate broadening parameter.  When using this, make sure that `log_gf` is the true
value (not adjusted for isotopic abundance).
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
For autoionizing lines (those for which E_upper > χ), Returns 0.0 for γ_vdW. Note the the "gamma"
values here are FWHM, not HWHM, of the Lorenztian component of the line profile.

In the calculation of `n*²`, uses the approximation that
``\\overbar{r^2} = 5/2 {n^*}^4 / Z^2``
which neglects the dependence on the angular momentum quantum number, l, in the the form given by
[Warner 1967](https://ui.adsabs.harvard.edu/abs/1967MNRAS.136..381W/abstract) (the earliest english
work reporting the Unsoeld result).
"""
function approximate_gammas(wl, species, E_lower; ionization_energies=ionization_energies)
    Z = species.charge + 1 #Z is ionization stage, not atomic number
    if ismolecule(species) || Z > 3
        return 0.0, 0.0
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

    Δrbar2 = (5 / 2) * Rydberg_eV^2 * Z^2 * (1 / (χ - E_upper)^2 - 1 / (χ - E_lower)^2)
    log_γvdW = if χ < E_upper
        0.0 # this will be interpretted as γ, rather than log γ, i.e. no vdW for autoionizing lines
    else
        # (log) γ_vdW From R J Rutten's course notes. 
        # Equations 11.29 and 11.30 from Gray 2005 are equivalent 
        6.33 + 0.4log10(Δrbar2) + 0.3log10(10_000) + log10(k)
    end

    γstark, log_γvdW
end

"""
    read_linelist(filename; format="vald", isotopic_abundances=Korg.isotopic_abundances)

Parse a linelist file, returning a vector of [`Line`](@ref)s.

The `format` keyword argument can be used to specify one of these linelist formats
(default: `"vald"`):

  - `"vald"` for a [VALD](http://vald.astro.uu.se/%7Evald/php/vald.php) linelist.
    These can be either "short" or "long" format,
    "extract all" or "extract stellar".  Air wavelengths will automatically be converted into vacuum
    wavelengths, and energy levels will be automatically converted from cm``^{-1}`` to eV.
  - `"kurucz"` for an atomic or molecular [Kurucz linelist](http://kurucz.harvard.edu/linelists.html)
    (format=kurucz_vac if it uses vacuum wavelengths; be warned that Korg will not assume that
    wavelengths are vacuum below 2000 Å),
  - `"moog"` for a [MOOG linelist](http://www.as.utexas.edu/%7Echris/moog.html)
    (doesn't support broadening parameters or dissociation energies, assumed to be in vacuum wavelengths).
  - `"moog_air"` for a MOOG linelist in air wavelengths.
  - `"turbospectrum"` for a
    [Turbospectrum linelist](https://github.com/bertrandplez/Turbospectrum2019/blob/master/DOC/Readme-Linelist_format_v.19)
    in air wavelengths. Note that Korg doesn't make use of the (optional) orbital angular momentum quantum number, l,
    for the upper or lower levels, so it won't fall back on generic ABO recipes when the ABO
    parameters are not available.
    Korg's interpretation of the `fdamp` parameter is also slightly different from Turbospectrum's.
    See the documentation of the `vdW` parameter of [`Line`](@ref) for details.  Korg will error if
    encounters an Unsoeld fudge factor, which it does not support.
  - "turbospectrum_vac" for a Turbospectrum linelist in vacuum wavelengths.

For VALD and Turbospectrum linelists with isotope information available, Korg will scale log gf
values by isotopic abundance (unless VALD has already pre-scaled them), using isotopic abundances
from [NIST](https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses)
([Korg.isotopic_abundances]).
To use custom isotopic abundances, just pass `isotopic_abundances` with the same structure:
a dict mapping atomic number to a dict mapping from atomic weight to abundance.

Be warned that for linelists which are pre-scaled for isotopic abundance, the estimation of
radiative broadening from log(gf) is not accurate.
"""
function read_linelist(fname::String; format="vald", isotopic_abundances=isotopic_abundances)
    format = lowercase(format)
    linelist = open(fname) do f
        if startswith(format, "kurucz")
            vac = endswith(format, "_vac")
            # open a new reader so we dont chomp the first line
            firstline = open(first_nonempty_line, fname)
            if length(firstline) > 100
                parse_kurucz_linelist(f; vacuum=vac)
            else
                parse_kurucz_molecular_linelist(f; vacuum=vac)
            end
        elseif format == "vald"
            parse_vald_linelist(f, isotopic_abundances)
        elseif format == "moog"
            parse_moog_linelist(f, isotopic_abundances, true)
        elseif format == "moog_air"
            parse_moog_linelist(f, isotopic_abundances, false)
        elseif format == "turbospectrum"
            parse_turbospectrum_linelist(f, isotopic_abundances, false)
        elseif format == "turbospectrum_vac"
            parse_turbospectrum_linelist(f, isotopic_abundances, true)
        else
            throw(ArgumentError("$(format) is not a supported linelist format"))
        end
    end

    filter!(linelist) do line #filter triply+ ionized and hydrogen lines
        (0 <= line.species.charge <= 2) && (line.species != species"H_I")
    end

    #ensure linelist is sorted
    if !issorted(linelist; by=l -> l.wl)
        sort!(linelist; by=l -> l.wl)
    end

    linelist
end

first_nonempty_line(io) =
    while !eof(io)
        line = readline(io; keep=true)
        if !all(isspace, line)
            return line
        end
    end

#used to handle missing gammas in vald and kurucz lineslist parsers
tentotheOrMissing(x) = x == 0 ? missing : 10^x
idOrMissing(x) = x == 0 ? missing : x

function parse_kurucz_linelist(f; vacuum=false)
    lines = Line{Float64,Float64,Float64,Float64,Float64,Float64}[]
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
            abs(parse(Float64, s)) * c_cgs * hplanck_eV
        end

        wl_transform = vacuum ? identity : air_to_vacuum

        push!(lines,
              Line(wl_transform(parse(Float64, row[1:11]) * 1e-7), #convert from nm to cm
                   parse(Float64, row[12:18]),
                   Species(row[19:24]),
                   min(E_levels...),
                   tentotheOrMissing(parse(Float64, row[81:86])),
                   tentotheOrMissing(parse(Float64, row[87:92])),
                   idOrMissing(parse(Float64, row[93:98]))))
    end
    lines
end

function parse_kurucz_molecular_linelist(f; vacuum=false)
    throw(ArgumentError("Kurucz linelists are not yet supported for molecules"))
    lines = Line[]
    for row in eachline(f)
        row == "" && continue #skip empty lines

        #kurucz provides wavenumbers for "level 1" and "level 2", which is which is 
        #determined by parity
        E_levels = map((row[23:32], row[39:48])) do s
            #abs because Kurucz multiplies predicted values by -1
            abs(parse(Float64, s)) * c_cgs * hplanck_eV
        end

        wl_transform = vacuum ? identity : air_to_vacuum

        push!(lines,
              Line(wl_transform(parse(Float64, row[1:10]) * 1e-7), #convert from nm to cm
                   parse(Float64, row[11:17]),
                   Species(row[49:52]),
                   min(E_levels...)))
    end
    lines
end

function parse_vald_linelist(f, isotopic_abundances)
    lines = filter!(collect(eachline(f))) do line
        length(line) > 0 && line[1] != '#' # remove comments and empty lines
    end

    # ignore truncation warning
    if startswith(lines[1], " WARNING: Output was truncated to 100000 lines")
        lines = lines[2:end]
    end

    lines = replace.(lines, "'" => "\"") #replace single quotes with double

    # is this an "extract all" or an "extract stellar" linelist?
    extractall = !occursin(r"^\s+\d", lines[1])
    firstline = extractall ? 3 : 4
    header = lines[firstline-1]

    scale_isotopes = any(startswith.(lines, "* oscillator strengths were NOT scaled "))
    if !scale_isotopes && !any(startswith.(lines, "* oscillator strengths were scaled "))
        throw(ArgumentError("Can't parse linelist.  Can't detect whether log(gf)s are scaled by " *
                            "isotopic abundance."))
    end

    #we take the linelist to be long-format when the second line after the header starts with a 
    #space or a quote followed a space.  In some linelists the quotes are there, but in others 
    #they are not.
    shortformat = !(occursin(r"^\"? ", lines[firstline+1]))
    body = lines[firstline:(shortformat ? 1 : 4):end]
    body = body[1:findfirst(l -> l[1] != '\"' || !isuppercase(l[2]), body)-1]

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
    body = CSV.File(reduce(vcat, codeunits.(body .* "\n")); header=CSVheader, delim=',',
                    silencewarnings=true)

    E_low = if contains(header, "cm") #convert E_low to eV if necessary
        body.E_low * c_cgs * hplanck_eV
    elseif contains(header, "eV")
        body.E_low
    else
        error("Can't parse linelist.  Can't determine energy units: " * E_col)
    end

    wl = 1e-8 * if contains(header, "air") #convert wls to vacuum if necessary
        air_to_vacuum.(body.wl)
    elseif contains(header, "vac")
        body.wl
    else
        error("Can't parse linelist.  Can't determine vac/air wls: " * header)
    end

    Δlog_gf = if scale_isotopes
        refs = if !shortformat #the references are on different lines
            # this line breaks the code formatter. 
            # https://github.com/domluna/JuliaFormatter.jl/issues/860
            #! format: off
            lines[firstline+3 .+ ((0:length(body)-1) .* 4)]
            #! format: on
        else #references are in the last column
            body.reference
        end

        map(refs) do ref
            #find things that look like (16)O or (64)Ni in reference string
            regexp = r"\((?<isotope>\d\d?\d?)\)(?<elem>\p{Lu}\p{Ll}?)"
            #add up the adjustments to log(gf) from isotopic abundances (zero if no info is present)
            log_probs = map(findall(regexp, ref)) do r
                m = match(regexp, ref[r])
                log10(isotopic_abundances[atomic_numbers[m["elem"]]][parse(Int64, m["isotope"])])
            end
            sum([0; log_probs])
        end
    else
        0
    end

    gamma_rad = map(wl, body.loggf, body.gamma_rad) do lambda, loggf, gamma
        if gamma == 0
            approximate_radiative_gamma(lambda, loggf)
        else
            10^gamma
        end
    end

    Line.(wl, body.loggf .+ Δlog_gf, Species.(body.species), E_low, gamma_rad,
          tentotheOrMissing.(body.gamma_stark),
          idOrMissing.(body.gamma_vdW))
end

#todo support moog linelists with broadening parameters?
function parse_moog_linelist(f, isotopic_abundances, vacuum_wavelengths)
    lines = collect(eachline(f))
    # The first line is ignored.  It's for human-readability only.
    linelist = map(lines[2:end]) do line
        toks = split(line)

        # special hanlding for the decimal part of species strings.  MOOG uses the first digit only
        # to represent the charge.  The rest contains isotopic info.
        dotind = findfirst('.', toks[2])
        spec = Species(toks[2][1:dotind+1])

        isostring = toks[2][dotind+2:end]
        #Note: this will fail if the atoms species code are not in order of atomic number.  This is always 
        #the case in the linelists I've seen.
        Δ_log_gf = if isostring == "" || !isnothing(match(r"^0+$", isostring)) #if all 0s
            0.0
        else
            natoms = length(get_atoms(spec))
            @assert length(isostring) % natoms == 0
            digits_per = length(isostring) ÷ natoms
            map(get_atoms(spec), 1:digits_per:length(isostring)-digits_per+1) do el, i
                m = parse(Int, isostring[i:i+digits_per-1])
                if m in keys(isotopic_abundances[el])
                    log10(isotopic_abundances[el][m])
                else
                    @info "No isotopic abundance for $(atomic_symbols[el]) $m. Leaving the log(gf) unchanged. (Occured when parsing $line)"
                    0.0
                end
            end |> sum
        end

        wl = parse(Float64, toks[1]) * 1e-8 #convert Å to cm
        if !vacuum_wavelengths
            wl = air_to_vacuum(wl)
        end

        Line(wl,
             parse(Float64, toks[4]) + Δ_log_gf,
             spec,
             parse(Float64, toks[3]))
    end
    linelist
end

function parse_turbospectrum_linelist(fn, isotopic_abundances, vacuum)
    # https://github.com/bertrandplez/Turbospectrum2019/blob/master/DOC/Readme-Linelist_format_v.19

    lines = readlines(fn)
    species_headers = filter(1:length(lines)) do i
        i != length(lines) && lines[i][1] == '\'' && lines[i+1][1] == '\''
    end

    transitions_for_each_species = map(1:length(species_headers)) do header_line_ind
        first_line_ind = species_headers[header_line_ind]
        last_line_ind = if header_line_ind == length(species_headers)
            length(lines)
        else
            species_headers[header_line_ind+1] - 1
        end

        # species line might look like this (carrot is beginning of line):
        # ^'  26.000            '    1       2342
        # here, the 26 refers to Fe (works as everything else does for molecules).  The decimal part
        # is the isotope information, NOT THE CHARGE.  The "1" is the ionization starge, i.e. the 
        # charge + 1. 2341 is the number of lines.

        species_line = lines[first_line_ind]
        m = match(r"'\s*(?<formula>\d+)\.(?<isostring>\d+)\s+'\s+(?<ion>\d+)\s+(?<n_lines>\d+)\s*",
                  species_line)
        formula = Formula(m["formula"])
        charge = parse(Int, m["ion"]) - 1
        spec = Korg.Species(formula, charge)
        n_lines = parse(Int, m["n_lines"])
        if last_line_ind - first_line_ind - 1 != n_lines
            error("Can't parse this line list.  The file says there are $n_lines lines for $spec, but I see $(last_line_ind - first_line_ind - 2) lines.")
        end

        isostring = m["isostring"]
        isotopic_Δ_loggf = map(get_atoms(spec), 1:3:length(isostring)-2) do el, i
            m = parse(Int, isostring[i:i+2])
            if m == 0 # no isotope specified for this constituent nucleus
                0.0
            else
                log10(isotopic_abundances[el][m])
            end
        end |> sum
        map(lines[first_line_ind+2:last_line_ind]) do line
            parse_turbospectrum_linelist_transition(spec, isotopic_Δ_loggf, line, vacuum)
        end
    end
    sort!(vcat(transitions_for_each_species...); by=l -> l.wl)
end

function parse_turbospectrum_linelist_transition(species, Δloggf, line, vacuum)
    # from the Turbospectrum docs (In practice linelists may have as few at 6 columns:
    #
    # For each line that follows:
    # col 1: lambda(A)  
    # col 2: Elow(eV) 
    # col 3: loggf 
    # col 4: fdamp (see below)
    # col 5: gup
    # col 6: gamma_rad (if =0, gf-value is used to compute gamma_rad)
    # col 7: gamma_Stark (may be omitted)
    # col 8: s,p,d,f etc for upper level (or X), see fdamp
    # col 9: same for lower level
    # col 10: equivalent width, when needed (abundance determination in eqwidt run)
    # col 11: error in eqw
    # col 12: (in quotes) some text describing levels or whatever you like to include

    # there could be a comma separating tokens (and fortran would parse), but I've never seen it.
    toks = split(line)

    log_gf = parse(Float64, toks[3])
    wl = air_to_vacuum(parse(Float64, toks[1]) * 1e-8)
    gamma_rad = parse(Float64, toks[6])
    if gamma_rad == 0 || gamma_rad == 1
        gamma_rad = Korg.approximate_radiative_gamma(wl, log_gf)
    end

    # if toks[7] is present, but gamma_stark is skipped, it will be the l for the upper level.
    stark_log_gamma = if length(toks) < 7 || isnothing(tryparse(Float64, toks[7]))
        missing
    else
        tentotheOrMissing(tryparse(Float64, toks[7]))
    end
    Elower = parse(Float64, toks[2])

    wltrans = vacuum ? identity : air_to_vacuum
    wl = wltrans(parse(Float64, toks[1]) * 1e-8)

    fdamp = parse(Float64, toks[4])

    Line(wl,
         log_gf + Δloggf,
         species,
         Elower,
         gamma_rad,
         stark_log_gamma,
         fdamp)
end

"""
    get_VALD_solar_linelist()

Get a VALD "extract stellar" linelist produced at solar parameters, with the "threshold" value
set to 0.01.  It was downloaded on 2021-05-20. It is intended to be used for quick tests only.
"""
get_VALD_solar_linelist() = read_linelist(joinpath(_data_dir, "linelists",
                                                   "vald_extract_stellar_solar_threshold001.vald"))

"""
    get_APOGEE_DR17_linelist(; include_water=true)

The APOGEE DR 17 linelist.  It ranges from roughly 15,000 Å to 17,000 Å.  It is
nearly the same at the DR 16 linelist described in
[Smith+ 2021](https://ui.adsabs.harvard.edu/abs/2021AJ....161..254S/abstract).
"""
function get_APOGEE_DR17_linelist(; include_water=true)
    dir = joinpath(_data_dir, "linelists", "APOGEE_DR17")

    atoms = read_linelist(joinpath(dir, "turbospec.20180901t20.atoms_no_ba");
                          format="turbospectrum")
    mols = read_linelist(joinpath(dir, "turbospec.20180901t20.molec"); format="turbospectrum")

    linelists = if include_water
        waterfile = joinpath(dir, "pokazatel_water_lines.h5")
        waterlines = Line.(Float64.(h5read(waterfile, "wl")),
                           Float64.(h5read(waterfile, "log_gf")),
                           Ref(species"H2O_I"),
                           Float64.(h5read(waterfile, "E_lower")),
                           Float64.(h5read(waterfile, "gamma_rad")))
        [atoms; mols; waterlines]
    else
        [atoms; mols]
    end

    sort!(linelists; by=l -> l.wl)
end

"""
    get_GALAH_DR3_linelist()

The GALAH [DR 3](https://www.galah-survey.org/dr3/overview/)
linelist (also used for [DR 4](https://github.com/svenbuder/GALAH_DR4)).
It ranges from roughly 4,675 Å to 7,930 Å.
This linelist is based on, but distinct from
[Heiter 2021](https://ui.adsabs.harvard.edu/abs/2021A%26A...645A.106H/).
See [Buder et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.506..150B%2F/abstract) for
details.
"""
function get_GALAH_DR3_linelist()
    path = joinpath(_data_dir, "linelists", "GALAH_DR3", "galah_dr3_linelist.h5")
    h5open(path, "r") do f
        species = map(eachcol(read(f["formula"])), read(f["ionization"])) do atoms, ion
            formula = if atoms[2] == 0x0
                Formula(atoms[1])
            elseif atoms[3] == 0x0
                Formula(sort(atoms[1:2]))
            else
                Formula(sort(atoms))
            end
            Species(formula, ion - 1)
        end
        lines = Line.(Float64.(read(f["wl"])),
                      Float64.(read(f["log_gf"])),
                      species,
                      Float64.(read(f["E_lo"])),
                      tentotheOrMissing.(Float64.(read(f["gamma_rad"]))),
                      tentotheOrMissing.(Float64.(read(f["gamma_stark"]))),
                      idOrMissing.(Float64.(read(f["vdW"]))))
        filter!(lines) do line
            #take out the hydrogen lines
            line.species != species"H I"
        end
    end
end

"""
    get_GES_linelist()

The Gaia-ESO survey linelist from
[Heiter et al. 2021](https://ui.adsabs.harvard.edu/abs/2021A&A...645A.106H/abstract).
"""
function get_GES_linelist()
    @warn "This function is may fail on some systems. See https://github.com/ajwheeler/Korg.jl/issues/309 for details."
    path = joinpath(artifact"Heiter_2021_GES_linelist",
                    "Heiter_et_al_2021_2022_06_17",
                    "Heiter_et_al_2021.h5")
    h5open(path, "r") do f
        Line.(Float64.(air_to_vacuum.(read(f["wl"]))),
              Float64.(read(f["log_gf"])),
              [Species(s) for s in read(f["species"])],
              Float64.(read(f["E_lower"])),
              tentotheOrMissing.(Float64.(read(f["gamma_rad"]))),
              tentotheOrMissing.(Float64.(read(f["gamma_stark"]))),
              read(f["vdW"]))
    end
end

"""
    _load_alpha_5000_linelist([path])

Load the default linelist for calculating the absorption coefficient at 5000 Å.  This for internal
use when the provided linelist doesn't cover the region and a radiative transfer scheme using
τ_5000 is used.

This linelist loaded into `Korg._alpha_5000_default_linelist` when Korg is imported.
"""
function _load_alpha_5000_linelist(path=joinpath(_data_dir, "linelists", "alpha_5000",
                                                 "alpha_5000_lines.csv"))
    csv = CSV.File(path)
    map(csv) do row
        vdW = if ',' in row.vdW
            σ, α = split(row.vdW[2:end-1], ',')
            (parse(Float64, σ), parse(Float64, α))
        else
            parse(Float64, row.vdW)
        end
        Line(row.wl, row.log_gf, Species(row.species), row.E_lower, row.gamma_rad, row.gamma_stark,
             vdW)
    end
end

"""
The default linelist for calculating the absorption coefficient at 5000 Å.  This for internal
use when the provided linelist doesn't cover the region and a radiative transfer scheme using
τ_5000 is used.

See also [`_load_alpha_5000_linelist`](@ref).
"""
const _alpha_5000_default_linelist = _load_alpha_5000_linelist()
