using Korg, HDF5, Base
using Statistics: mean
using Interpolations

include("NIST_NORAD_TOPBase_parsing_code.jl")

"""
Returns cross section in megabarns summed from all electron configurations for a given species.

TODO
"""
function single_species_bf_cross_section(spec::Species, λs, Ts, data_dir)
    #this method handles loading the correct files if you just specify the species
    @assert !ismolecule(spec)

    Z = get_atoms(spec.formula)[1]
    n_electrons = Z - spec.charge

    filename = "p" * lpad(Z, 2, "0") * "." * lpad(n_electrons, 2, "0") * ".dat"
    #use NORAD tables when available, TOPBase otherwise
    path, norad_format = if isfile(joinpath(data_dir, "bf_cross_sections", "NORAD", filename))
        joinpath(data_dir, "bf_cross_sections", "NORAD", filename), true
    else
        joinpath(data_dir, "bf_cross_sections", "TOPBase", filename), false
    end
    cross_sections = parse_TOPBase_cross_sections(path, norad_format)

    nist_levels = parse_NIST_energy_levels(joinpath(data_dir, "NIST_energy_levels"), spec)

    Z = get_atoms(spec)[1]
    ionization_energy = Korg.ionization_energies[Z][spec.charge + 1]

    single_species_bf_cross_section(cross_sections, nist_levels, ionization_energy, 
                             Korg.partition_funcs[spec], λs, Ts)
end
function single_species_bf_cross_section(cross_sections, energy_levels, ionization_energy, U, λs, Ts)
    # convert λ_vals to photon energies in eV
    photon_energies = (hplanck_eV * c_cgs) ./ λs

    # precompute Temperature-dependent constant
    Us = U.(log.(Ts))
    β = 1.0 ./ (kboltz_eV .* Ts)
    
    #ionization_energy_topbase = maximum(E for (E, ) in values(cross_sections))
    #display("ionization energy: $(ionization_energy) (NIST) $(ionization_energy_topbase) (theoretical)")

    # prepare the output array where results will be accumulated.  This the cross section obtained 
    # by summing over electron states
    total_sigma = zeros(eltype(photon_energies), (length(photon_energies),length(Ts)))
    total_weight = zeros(length(Ts))

    for (state, (topbase_binding_energy, Es, σs)) in pairs(cross_sections)
        #the topbase binding energies are not precise.  Use NIST empirical energy levels instead.
        energy_level = if state in keys(energy_levels)
            energy_levels[state]
        else
            #println("skipping $state")
            continue
        end

        empirical_binding_energy = ionization_energy - energy_level #in eV
        #display("binding energy: $(empirical_binding_energy) (NIST) $(topbase_binding_energy) (theoretical)")

        # adjust Es (photon energies) to match empirical binding energy
        Es .+= empirical_binding_energy - topbase_binding_energy
        Interpolations.deduplicate_knots!(Es, move_knots=true) #shift repeat Es to the next float

        σ_itp = LinearInterpolation(Es, σs, extrapolation_bc=0.0)

        #ion_energies are the energies with respect to the ionization energy of the species
        #excitation_potential_ryd = ionization_energy - binding_energy

        # g*exp(-βε)/U at each temperature
        weights = statistical_weight(state) .* exp.(-energy_level .* β) ./ Us
        display(weights)
        total_weight += weights

        # stimulated emission is rare since it requires a free electron with the right energy
        # the extra factor of nₑ makes it negligable?
        #total_sigma .+= σ_itp.(photon_energies) .* (1.0 .- exp.(-photon_energies.*β_Ryd')) .* weights'
        total_sigma .+= σ_itp.(photon_energies) .* weights'
    end

    display(total_weight)

    total_sigma 
end

function save_bf_cross_section(f, spec, logTs, νs, ξ=2e5)
    @assert issorted(νs, rev=true) # should be decreasing freq
    
    Ts = 10 .^ logTs
    wls = Korg.c_cgs ./ νs
    σs = single_species_bf_cross_section(spec, wls, Ts, DATA_DIR)
    
    Rs = Korg.c_cgs ./ sqrt.(2 * Korg.kboltz_cgs * Ts ./ Korg.get_mass(Korg.species"He I") .+ ξ^2) 
    λs = Korg.c_cgs ./ νs # the frequencecy grid is linear-enough in wavelength over the ~ 1 Å scale LSF 

    # cast to single precision
    σs = hcat(Korg.constant_R_LSF.(eachcol(σs), Ref(λs), Rs)...)
    mask = σs .< floatmin(Float32)
    σs[mask] .= 0
    print(spec)
    println(" ", minimum(σs), " -- ", maximum(σs), " ($(mean(mask)) denormal)")
    log_σs = Float32.(log.((σs)))
    
    # save in order of increasing freq / decreasing wl
    f["cross-sections/"*string(spec), shuffle=(), deflate=6] = reverse(log_σs, dims=1)
end

DATA_DIR = ARGS[1]
OUT_FILE = "bf_cross-sections.h5"

println("writing $OUT_FILE")

logTs = 2:0.1:5
wl_lo = 490 * 1e-8
wl_hi = 30_050 * 1e-8

# descending freq, ascending wl
νs = Korg.c_cgs / wl_lo : -1e11 : Korg.c_cgs / wl_hi


filename = OUT_FILE
h5open(filename, "w") do f
    f["logT_min"]  = logTs[1]
    f["logT_step"] = step(logTs)
    f["logT_max"] = logTs[end]
    f["nu_min"] = νs[end]
    f["nu_step"] = -step(νs)
    f["nu_max"] = νs[1]

    atomic_numbers = [6, 12, 13, 14]
    #atomic_numbers = [1:14... ; 16:2:20... ; 26]
    @time for Z in atomic_numbers
        #for ionization_number in [1, 2]
        for ionization_number in [1]
            if Z == 1 && ionization_number == 2
                continue # no such thing as H II 
            end
            spec = Korg.Species(Korg.atomic_symbols[Z] * " $ionization_number")
            save_bf_cross_section(f, spec, logTs, νs)
        end
    end
end
;
