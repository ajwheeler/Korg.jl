using Base
using Interpolations: LinearInterpolation, deduplicate_knots!
using Korg: partition_funcs, Species, ismolecule, get_atoms, _data_dir, move_bounds, @species_str
using Korg: hplanck_eV, c_cgs, RydbergH_eV, kboltz_eV # constants

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

function Base.show(io::IO, state::ElectronState)
    orbital_letter = if state.L == 0
        "S"
    elseif state.L == 1
        "P"
    else
        "DFGHIJKLMNO"[state.L - 1] #2->D, 3->F, etc.
    end
    parity_char = if state.P == 0
        ""
    else
        "*"
    end
    print(io, state.spin_multiplicity, orbital_letter, parity_char, " (level ", state.iLV, ")")
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

function parse_TOPBase_cross_sections(filename, norad_format=false)
    cross_sections = Dict{ElectronState, Any}()
    lines = readlines(filename)

    #points to "current" line. This could be trivially rewritten to stream line-by-line from disk.
    i = 1 

    #skip the big header in norad formatted files
    if norad_format
        while (lines[i] != "-----------------------------------------------------------------------")
            i += 1
        end
        i += 1
    end

    i += 1 # the first line tells you which species the file corresponds to, skip to the second

    while (i <= length(lines)) && (!startswith(lines[i], _FINAL_LINE))
        #parse the header specifying the state
        state = ElectronState(parse(UInt8, lines[i][4:5]), parse(UInt8, lines[i][9:10]),
                              parse(UInt8, lines[i][14:15]), parse(UInt8, lines[i][19:20]))

        #the next line tells you how many points there are
        # for example: "     479    489" indicates 489 points. (The "479" is an internal thing.)
        npoints = parse(Int, lines[i+1][9:14])
        # The line after that indicates the binding energy.
        # for example: "  1.000000E+00    0.0100" indicates that the binding energy is 1 Ryd.
        binding_energy = parse(Float32, lines[i+2][3:14])

        # map the electron state the cross-section interpolator in the cross_sections dict.
        Es = Vector{Float64}(undef, npoints)
        σs = Vector{Float64}(undef, npoints)

        i += 3 # Make i point at the first real line of the cross-sections vals
        for j in 1:npoints
            Es[j] = parse(Float32, lines[i+j-1][3:14])
            σs[j] = parse(Float32, lines[i+j-1][16:24])
        end
        deduplicate_knots!(Es, move_knots=true) # shift repeated E values to the next float 
        itp = LinearInterpolation(Es, σs, extrapolation_bc=0.0)

        cross_sections[state] = binding_energy, itp
        i += npoints #move i to the next electron state
    end

    cross_sections
end

"""
TODO
"""
function cross_section_bf(cross_sections, U, λs, Ts)
    # convert λ_vals to photon energies in Ryd
    photon_energies = (hplanck_eV * c_cgs / RydbergH_eV) ./ λs

    # precompute Temperature-dependent constant
    Us = U.(log.(Ts))
    β_Ryd = RydbergH_eV./(kboltz_eV .* Ts)
    binding_energy_ground = maximum(E for (E, itp) in values(cross_sections))

    # prepare the output array where results will be accumulated.  This the cross section obtained 
    # by averaging over electron states
    weighted_average = zeros(eltype(photon_energies), (length(photon_energies),length(Ts)))

    for (state, (binding_energy, cross_section_itp)) in pairs(cross_sections)
        #ion_energies are the energies with respect to the ionization energy of the species
        excitation_potential_ryd = binding_energy_ground - binding_energy

        # g*exp(-βε)/U at each temperature
        weights = statistical_weight(state) .* exp.(-excitation_potential_ryd.*β_Ryd) ./ Us

        #cross section for this state including stimulated emission term at each λ
        σs = cross_section_itp.(photon_energies) .* (1.0 .- exp.(-photon_energies.*β_Ryd'))

        weighted_average .+= σs .* weights'
    end

    weighted_average * 1e-18 #megabarns -> cm^2
end
function cross_section_bf(spec::Species, λs, Ts, data_dir)
    @assert !ismolecule(spec)

    #use NORAD tables for iron, TOPBase for others
    use_norad = (spec == species"Fe I") || (spec == species"Fe II")
    dir = joinpath(data_dir, use_norad ? "NORAD" : "TOPBase")

    Z = get_atoms(spec.formula)[1]
    n_electrons = Z - spec.charge
    filename = joinpath(dir, "p" * lpad(Z, 2, "0") * "." * lpad(n_electrons, 2, "0") * ".dat")
    cross_sections = parse_TOPBase_cross_sections(filename, use_norad)

    cross_section_bf(cross_sections, partition_funcs[spec], λs, Ts)
end
