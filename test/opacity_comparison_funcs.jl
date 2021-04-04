# This file defines some functions that are used to compare of our some opacity calculations
# against the values shown in several panels of Figure 8.5 from Gray's 2005 edition of "The
# Observation and Analysis of Stellar Photospheres".
#
# Data was extracted using the WebPlotDigitizer tool. We've recorded values with greater
# precision than the values should actually have. The values are primarily meant to be used to
# make sure that our implementation give results to the correct order of magniutde and
# qualitatively has the expected wavelength dependences.
#
# Tests are explicitly NOT defined in this file so that functions defined in this file can be used
# in scripts outside of this framework (to create plots in order to assess disagreement)

using HDF5
using Interpolations: LinearInterpolation, Throw

# Load Gray05 data for a given panel. This returns a tuple holding two dictionaries:
# 1. The first holds the panel properties (i.e. the temperature and partial electron pressure that
#    was assumed for the calculation)
# 2. The second holds a 2-tuple for each source of opacity source. The first tuple entry specifies
#    wavelengths (in Ångstroms) and the second entry holds the opacity contribution at those
#    wavelengths.
function load_panel_data(panel_name, h5fname = joinpath(@__DIR__, "data/gray05_fig8.5.h5"))
    attribute_dict = Dict{String,Float64}()
    data_dict = Dict{String,Tuple{Array{Float64,1},Array{Float64,1}}}()

    h5open(h5fname, "r") do file
        # at the highest level, the file is divided into one group per panel
        panel_group = file[panel_name]

        # load the attributes from the panel group. This nominally holds the temperature and
        # partial electron pressure for which the values are computed.
        for key in keys(attributes(panel_group))
            attribute_dict[key] = read(attributes(panel_group)[key])
        end

        # In a given panel_group there is one subgroup per opacity source.
        for opacity_subgroup_name in keys(panel_group)
            opacity_subgroup = panel_group[opacity_subgroup_name]
            if !(length(keys(opacity_subgroup)) == 2 &&
                 haskey(opacity_subgroup, "lambda") &&
                 haskey(opacity_subgroup, "data"))
                msg = string("Each opacity subgroup is expected to just have a \"lambda\" ",
                             "dataset and a \"data\" dataset") 
                throw(ErrorException(msg))
            end
            tmp = (read(opacity_subgroup["lambda"]), read(opacity_subgroup["data"]))
            if length(tmp[1]) != length(tmp[2])
                throw(ErrorException("Datasets within a given subgroup should have a fixed length"))
            end
            data_dict[opacity_subgroup_name] = tmp
        end
    end
    (attribute_dict, data_dict)
end

# To facillitate comparisons against opacity from H₂⁺ and He⁻, we need to have reasonable estimates
# of the ratio of free electrons to protons. In the future, we might want to hardcode the
# calculated values so that changes to the atomic data doesn't break all of these tests.

# from table D.2. This interpolates the log of the partition function.
# I updated cases were the table records 0.301 with log10(2.0)
function _gray05_H_I_partition_func(T)
    log_2 = log10(2.0)
    θ_vals = [  0.2,   0.4,   0.6,   0.8,   1.0,   1.2,   1.4,   1.6,   1.8,   2.0]
    table  = [0.368, 0.303, log_2, log_2, log_2, log_2, log_2, log_2, log_2, log_2]
    interpolator = LinearInterpolation(θ_vals, table, extrapolation_bc=Throw())
    10^interpolator(5040.0 / T)
end

const gray05_partition_funcs =
    Dict( "H_I" => _gray05_H_I_partition_func,
          "H_II" => T -> 1.0,
          "He_I" => T -> (@assert 2520 ≤ T ≤ 25200, 1.0),
          "He_II" => T -> (@assert 2520 ≤ T ≤ 25200, 2.0),
          "He_III" => T -> 1.0)

# compute the ratio of the nₑ to n_H (the number density of all H I and H II)
# this function could be reused in the future for testing solutions to coupled Saha equations
function free_electrons_per_Hydrogen_particle(nₑ, T, abundances = SSSynth.solar_abundances,
                                              ionization_energies = SSSynth.ionization_energies,
                                              partition_funcs = SSSynth.partition_funcs)
    out = 0.0
    for element in SSSynth.atomic_symbols

        (χs, Us) = if element == "H"
            (ionization_energies[element][1:2], [partition_funcs["H_I"], partition_funcs["H_II"]])
        elseif haskey(ionization_energies,element)
            (ionization_energies[element],
             [partition_funcs[string(element, "_I")],
              partition_funcs[string(element, "_II")],
              partition_funcs[string(element, "_III")]])
        else
            continue
        end

        ion_state_weights = SSSynth.saha(χs, Us, T, nₑ)
        nₑ_per_ndens_species = 0.0
        for (i,ion_state_weight) in enumerate(ion_state_weights)
            num_electrons_from_state = i - 1
            nₑ_per_ndens_species += num_electrons_from_state * ion_state_weight
        end

        abundance = 10.0^(abundances[element]-12.0)
        out += abundance * nₑ_per_ndens_species
    end
    out
end

# At a given Temperature and number density of electrons, generate semi-realistic values for:
#    - number densities of H I
#    - number density of all hydgrogen particles (H I and H II)
#    - the mass density
#
# This is mostly used for comparisons against panels of figure 8.5 from Gray (2005). For most
# curves, we can simply assert that there is a fixed ionization fraction because we divide out the
# dependence on all of these terms during the calculation
function _semi_realisitic_dens(ne::Flt, fion::Flt = 0.02,
                               HydrogenMassFrac::Flt = 0.76) where {Flt<:AbstractFloat}
    # In the interest of semi-realistic numbers, Let's assume all free electrons come from Hydrogen
    nH_I = ne/fion
    nH = ne + nH_I
    ρ = nH * 1.67e-24/HydrogenMassFrac
    (nH_I, nH, ρ)
end

# Now actually define the functions that compute the opacities in the form comparable with Gray05

# Combined H I bound-free and free-free opacity
function HI_coefficient(λ, T, Pₑ, H_I_ion_energy = 13.598)
    # λ should have units of Ångstroms
    # Pₑ is the partial pressure of the electrons in dyne/cm²
    nₑ =  Pₑ/(SSSynth.kboltz_cgs * T)
    ν = (SSSynth.c_cgs*1e8)/λ

    bf_coef = begin
        H_I_partition_val = 2.0 # implicitly in the implementation provided by Gray (2005)
        nH_I = nₑ * 100.0 # this is totally arbitrary
        ρ = nH_I * 1.67e-24/0.76 # this is totally arbitrary
        bf_opac = SSSynth.ContinuumOpacity.H_I_bf(nH_I/H_I_partition_val, ν, ρ, T, H_I_ion_energy)
        bf_opac * ρ / (Pₑ * nH_I)
    end

    ff_coef = begin
        χs = [H_I_ion_energy, -1.0]
        Us = [T -> 2.0, T -> 1.0] # assumption in the approach described by Gray (2005)
        nH_total = nₑ * 100.0 # this is totally arbitrary
        ρ = nH_total * 1.67e-24/0.76 # this is totally arbitrary

        weights = SSSynth.saha(χs, Us, T, nₑ)
        nH_I = nH_total * weights[1]
        nH_II = nH_total * weights[2]
        ff_opac = SSSynth.ContinuumOpacity.H_I_ff(nH_II, nₑ, ν, ρ, T)
        ff_opac * ρ / (Pₑ * nH_I)
    end

    bf_coef + ff_coef
end

function Hminus_bf_coefficient(λ, T, Pₑ, ion_energy_H⁻ = 0.7552)
    # λ should have units of Ångstroms
    # Pₑ is the partial pressure of the electrons in dyne/cm²
    ne =  Pₑ/(SSSynth.kboltz_cgs * T)

    # the values of ndens, ρ and nH_I shouldn't actually matter since we will
    # just divide out all dependence on these variables.
    nH_I, nH, ρ = _semi_realisitic_dens(ne)

    partition_func = 2.0 # may want to include the temperature dependence of the partition function
    ν = (SSSynth.c_cgs*1e8)/λ

    opacity = SSSynth.ContinuumOpacity.Hminus_bf(nH_I/partition_func, ne, ν, ρ, T, ion_energy_H⁻)
    opacity * ρ / (Pₑ * nH_I)
end

function Hminus_ff_coefficient(λ, T, Pₑ)
    # λ should have units of Ångstroms
    # Pₑ is the partial pressure of the electrons in dyne/cm²
    ne =  Pₑ/(SSSynth.kboltz_cgs * T)

    # the values of ndens, ρ and nH_I shouldn't actually matter since we will
    # just divide out all dependence on these variables.
    nH_I, nH, ρ = _semi_realisitic_dens(ne)

    partition_func = 2.0 # may want to include the temperature dependence of the partition function
    ν = (SSSynth.c_cgs*1e8)/λ

    opacity = SSSynth.ContinuumOpacity.Hminus_ff(nH_I/partition_func, ne, ν, ρ, T)
    opacity * ρ / (Pₑ * nH_I)
end

# computes the combine H₂⁺ free-free and bound-free absorption in units of cm^2 per H atom (not a
# typo)
# There seems to be an error in our function to compute H₂⁺
function H2plus_coefficient(λ, T, Pₑ)
    # λ should have units of Ångstroms
    # Pₑ is the partial pressure of the electrons in dyne/cm²
    ne =  Pₑ/(SSSynth.kboltz_cgs * T)

    # in the future, we may want to use the abundances, ionization_energies and partition functions
    # used in Gray (2005) since we are comparing to his plots
    ne_div_nH = 
        free_electrons_per_Hydrogen_particle(ne, T)
    nH = ne/ne_div_nH

    χs = SSSynth.ionization_energies["H"][1:2]
    Us = [SSSynth.partition_funcs["H_I"], SSSynth.partition_funcs["H_II"]]
    weights = SSSynth.saha(χs, Us, T, ne)
    nH_I = weights[1]*nH
    nH_II = weights[2]*nH

    # set the partition function to 2.0 in nH_I_div_partition so that n(H I, n=1) = n(H I) for this
    # calculation, which is an assumption that Gray (2005) implicitly uses
    nH_I_div_partition = nH_I/2.0
    ρ = 1.0 # arbitrary value because we divide it out after

    ν = (SSSynth.c_cgs*1e8)/λ
    opacity = SSSynth.ContinuumOpacity.H2plus_bf_and_ff(nH_I_div_partition, nH_II, ν, ρ, T)
    opacity * ρ / (Pₑ * nH_I)
end



# compute He⁻ free-free absorption in units of cm^2 per H atom (not a typo)
function Heminus_ff_coefficient(λ, T, Pₑ)
    # λ should have units of Ångstroms
    # Pₑ is the partial pressure of the electrons in dyne/cm²
    ne =  Pₑ/(SSSynth.kboltz_cgs * T)

    # the values of ndens, ρ and nH_I shouldn't actually matter since we will
    # just divide out all dependence on these variables.

    # In the interest of semi-realistic numbers, Let's assume all free
    # electrons come from Hydrogen
    nH_I, nH, ρ = _semi_realisitic_dens(ne)

    # we need to compute the number density of He I. Since this function exists to facillitate
    # comparisons against plots from Gray (2005), we should employ the abundances from
    # their book. They record that there is 8.51e-2 He particles per H particle in table 16.3
    nHe = nH*8.51e-2
    # we need to figure out the abundance of neutral Helium. We will use the ionization potentials
    # tabulated in table D.1
    χs = [24.587]
    # technically, this should be χs = [24.587, 54.418] to account for both single and double
    # ionization, but equation 8.16 in Gray (2005) only accounts for single ionization.
    Us = [x -> 1.0,      # Table D.2 of Gray (2005) records the partion function of He II to be
          x -> 10^0.301] # as 10^0.301, the closest value at the table's precision to 2.0
                         
    weights = SSSynth.saha(χs, Us, T, ne)
    nHe_I = weights[1]*nHe

    # may want to include the temperature dependence of the partition function
    partition_func = 2.0
    ν = (SSSynth.c_cgs*1e8)/λ


    opacity = SSSynth.ContinuumOpacity.Heminus_ff(nHe_I/Us[1](T), ne, ν, ρ, T)
    opacity * ρ / (Pₑ * nH_I)
end

struct Bounds
    lower::Union{Float64, Nothing}
    upper::Union{Float64, Nothing}
    function Bounds(lower::Union{Float64, Nothing}, upper::Union{Float64, Nothing})
        if (!isnothing(lower)) && (!isnothing(upper)) && (lower >= upper)
            error("lower exceeds upper")
        end
        new(lower,upper)
    end
end

function inbounds(bounds::Bounds, vals::Array{T}) where {T<:AbstractFloat}
    map(vals) do val
        ((isnothing(bounds.lower) || bounds.lower <= val) &&
         (isnothing(bounds.upper) || bounds.upper >= val))
    end |> Array{Bool}
end

# There appears to be some errors in the H₂⁺ opacities, skipping them for now
const Gray05_opacity_form_funcs =
    Dict("H"          => (HI_coefficient,         Bounds(nothing, nothing), "H I bf and ff"),
         "Hminus_bf"  => (Hminus_bf_coefficient,  Bounds(2250.0, 15000.0),  "H⁻ bound-free"),
         "Hminus_ff"  => (Hminus_ff_coefficient,  Bounds(2604.0, 113918.0), "H⁻ free-free"),
         "Heminus_ff" => (Heminus_ff_coefficient, Bounds(5063.0, 151878.0), "He⁻ free-free"),
         "H2plus"     => (H2plus_coefficient,     Bounds(3847.0, 25000.0),   "H₂⁺ ff and bf"),
         )

# the absolute tolerances for values the opacity contributions from
#    - H⁻ bound-free
#    - H⁻ free-free
#    - He⁻ free-free
# There is a much larger discrepancy between the summed H I bf and ff opacity at large λ. But
# unlike for H₂⁺, this could simply be a consequence of the different approximaitons that are used
const Gray05_atols = Dict("a" => 0.05, "b" => 0.02, "c" => 0.04)

# this function is defined to make tests easier
function Gray05_comparison_vals(panel, opacity_func_name)
    panel_prop, panel_dict = load_panel_data(panel)
    orig_λ_vals, ref_data = panel_dict[opacity_func_name]
    temperature = panel_prop["temperature"]
    Pₑ = 10^panel_prop["log10(electron_pressure)"]

    func, bounds = Gray05_opacity_form_funcs[opacity_func_name][1:2]

    w = inbounds(bounds,orig_λ_vals)
    calculated_vals = func.(orig_λ_vals[w], temperature, Pₑ)
    (calculated_vals*1e26, ref_data[w])
end
