using Base
using Interpolations

"""
TODO
"""
struct ElectronState
    spin_multiplicity::UInt8 # 2S + 1
    L::UInt8                 # total orbital angular momentum 
    P::UInt8                 # parity 0->even, 1->odd
    iLV::UInt8               # index of the level within its spectroscopic series 
    function ElectronState(spin_multiplicity, L, P, iLV)
        @assert 1 <= spin_multiplicity <= 8 #TODO explain
        @assert P == 0 || P == 1            # parity must be either even or odd
        @assert iLV > 0                     # these are 1-indexed
        new(spin_multiplicity, L, P, iLV)
    end
end

"""
TODO
"""
function statistical_weight(state::ElectronState)
    state.spin_multiplicity * (2*state.L + 1)
end

_FINAL_LINE =  "    0    0    0    0"

# TODO
# In a couple tables, there's a place where the photon energy is listed out of order. This seems to
# be okay given that this is in a section of the table where the cross-section is zero.
const _species_with_problematic_subtables = (Species("He_I"), Species("C_I"), Species("C_II"),
                                             Species("Al_II"), Species("O_II"))


function parse_TOPBase_energy_levels(filename)
    energies = Dict{ElectronState, Float32}()
    lines = readlines(filename)

    i = 2 # the first line tells you which species the file corresponds to, skip to the second
    while (i <= length(lines)) && (!startswith(lines[i], _FINAL_LINE))
        #file is dived into sections by spectroscopic series, with each headed by a line like this:
        # "    3    2    0    8" 
        #grab each number with a regular expression
        series_match = match(r" +(?<spin_multiplicity>\d\d?) +(?<L>\d\d?) +(?<parity>\d\d?) +(?<n_levels>\d\d?)", lines[i])
        n_levels = parse(Int, series_match["n_levels"])

        #we ignore the nexline, the ones after that each correspond to an energy level of the series
        for j in i+2 : i+n_levels+1
            # These lines look like this: "    1  T  1  3 2  -1.11272E-0"
            # the first number is iLV and the last is the energy level

            #now map the state to the energy level the energies dict
            state = ElectronState(
                parse(UInt8, series_match["spin_multiplicity"]),
                parse(UInt8, series_match["L"]),
                parse(UInt8, series_match["parity"]),
                parse(UInt8, lines[j][4:5]))
            E = parse(Float32, lines[j][19:30])
            energies[state] = E
        end
        i += n_levels+2 # jump to next series
    end
    energies
end

function parse_TOPBase_cross_sections(filename)
    cross_sections = Dict{ElectronState, Any}()
    lines = readlines(filename)

    i = 2 # the first line tells you which species the file corresponds to, skip to the second
    while (i <= length(lines)) && (!startswith(lines[i], _FINAL_LINE))
        #parse the header specifying the state
        state = ElectronState(parse(UInt8, lines[i][4:5]), parse(UInt8, lines[i][9:10]),
                              parse(UInt8, lines[i][14:15]), parse(UInt8, lines[i][19:20]))

        #the next line tells you how many points there are
        # for example: "     479    489" indicates 489 points
        npoints = parse(Int, lines[i+1][9:14])

        # ignore the next line, it seems to exist to force the interpolator down to 0 and it doesn't 
        # count towards npoints. Make i point at the first real line of the cross-sections vals
        i += 3 

        # map the electron state the cross-section interpolator in the cross_sections dict.
        Es = Vector{Float32}(undef, npoints)
        σs = Vector{Float32}(undef, npoints)
        for j in 1:npoints
            Es[j] = parse(Float32, lines[i+j-1][3:14]) 
            σs[j] = parse(Float32, lines[i+j-1][16:24])
        end
        cross_sections[state] = CubicSpline(Es, σs)

        i += npoints #move i to the next electron state
    end

    cross_sections
end

"""
TODO
"""
function cross_section_bf_TOPBase(σfile, Efile,  U, χ,  λs, Ts)
    ion_energies = parse_TOPBase_energy_levels(Efile)
    cross_sections = parse_TOPBase_cross_sections(σfile)

    @assert keys(ion_energies) == keys(cross_sections)

    # convert λ_vals to photon energies
    photon_energies = (hplanck_eV * c_cgs / RydbergH_eV) ./ λs

    # precompute Temperature-dependent constant
    Us = U.(log.(Ts))
    β_Ryd = RydbergH_eV./(kboltz_eV .* Ts)

    # prepare the output array where results will be accumulated.  This the cross section obtained 
    # by averaging over electron states
    weighted_average = zeros(eltype(photon_energies), (length(photon_energies),length(Ts)))

    for state in keys(ion_energies)
        #ion_energies are teh energies with respect to the ionization energy of the species
        excitation_potential_ryd = χ + ion_energies[state] 

        # g*exp(-βε)/U at each temperature
        weights = statistical_weight(state) .* exp.(-excitation_potential_ryd.*β_Ryd) ./ Us
        #println(state)
        #println(ion_energies[state])
        #println(excitation_potential_ryd)
        #println(weights)

        #cross section for this state including stimulated emission term at each λ
        σs = cross_sections[state].(photon_energies; oob_val=0.0) .* (1.0 .- exp.(-photon_energies.*β_Ryd'))

        weighted_average .+= σs .* weights'
    end
    weighted_average * 1e-18 #megabarns -> cm^2
end
function cross_section_bf_TOPBase(spec::Species, λs, Ts; dir=joinpath(_data_dir, "TOPBaseData"))
    @assert !ismolecule(spec)
    Z = get_atoms(spec.formula)[1]
    n_electrons = Z - spec.charge
    file_suffix = lpad(Z, 2, "0") * "." * lpad(n_electrons, 2, "0") * ".dat"
    χ = ionization_energies[Z][spec.charge+1] / Rydberg_eV
    cross_section_bf_TOPBase(joinpath(dir, "p"*file_suffix), joinpath(dir, "e"*file_suffix),
                             partition_funcs[spec], χ, λs, Ts)
end
