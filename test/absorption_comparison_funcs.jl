# This file defines different sets of functions that are used for comparisons of some of our
# opacity calculations. Each set of functions is organized into it's own submodule.
#
# Tests are explicitly NOT defined in this file so that functions defined in this file can be used
# in scripts outside of this framework (to create plots in order to assess disagreement)

# ================================================================================================

# The first set of functions are defined to facillitate comparisons against the values shown in
# several panels of Figure 8.5 from Gray's 2005 edition of "The Observation and Analysis of Stellar
# Photospheres".

# For these panels, data was extracted using the WebPlotDigitizer tool. We've recorded values with
# greater precision than the values should actually have. The values are primarily meant to be used
# to make sure that our implementation give results to the correct order of magniutde and
# qualitatively has the expected wavelength dependences.

module Gray_opac_compare
using Korg, HDF5
using Korg: Interval, closed_interval, contained
using Interpolations: linear_interpolation, Throw

# Load Gray05 data for a given panel. This returns a tuple holding two dictionaries:
# 1. The first holds the panel properties (i.e. the temperature and partial electron pressure that
#    was assumed for the calculation)
# 2. The second holds a 2-tuple for each source of opacity source. The first tuple entry specifies
#    wavelengths (in Ångstroms) and the second entry holds the opacity contribution at those
#    wavelengths.
function load_panel_data(panel_name, h5fname=joinpath(@__DIR__, "data/gray05_fig8.5.h5"))
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
    θ_vals = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    table = [0.368, 0.303, log_2, log_2, log_2, log_2, log_2, log_2, log_2, log_2]
    interpolator = linear_interpolation(θ_vals, table; extrapolation_bc=Throw())
    10^interpolator(5040.0 / T)
end

# compute the ratio of the nₑ to n_H (the number density of all H I and H II)
# this function could be reused in the future for testing solutions to coupled Saha equations
function free_electrons_per_Hydrogen_particle(nₑ, T, abundances=Korg.default_solar_abundances)
    out = 0.0
    for element in 1:Korg.MAX_ATOMIC_NUMBER
        wII, wIII = Korg.saha_ion_weights(T, nₑ, element, Korg.ionization_energies,
                                          Korg.default_partition_funcs)

        nₑ_per_ndens_species = wII / (1 + wII + wIII) + 2wIII / (1 + wII + wIII)
        abundance = 10.0^(abundances[element] - 12.0)
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
function _semi_realisitic_dens(ne::F, fion::F=0.02, HydrogenMassFrac::F=0.76) where F<:Real
    # In the interest of semi-realistic numbers, Let's assume all free electrons come from Hydrogen
    nH_I = ne / fion
    nH = ne + nH_I
    ρ = nH * 1.67e-24 / HydrogenMassFrac
    (nH_I, nH, ρ)
end

# Gray05 uses some confusing terminology and notation.
# - he defines the "atomic absorption coefficient" as the area per absorber and represent it with
#   the variable α. For reference, Rybicki & Lightman call this same quantity the cross-section and
#   represent it with the variable σ_ν.
# - he defines the "continuous absorption coefficient per neutral hydrogen atom" and denotes it
#   with κ. Rybicki & Lightman don't explicitly define this quantity, but it's compatible with
#   their nomenclature. R&L would express this quantity as α_ν/nH_I (in their notation). This is
#   NOT a cross-section in the general case (e.g. for continuum absorption of He, Gray still
#   defines the "continuous absorption coefficient per neutral hydrogen atom" as κ = α_ν/nH_I).
# - he uses κ_ν to denote the mass opacity (just like Rybicki & Lightman)
#
# Now we actually define the functions that compute the opacities in the form comparable with
# Gray05. These functons return the "continuous absorption coefficient per neutral hydrogen atom"
# divided by Pₑ, the partial electron pressure.
#
# For the sake of clarity, we label Rybicki & Lightman's α_ν = κ_ν*ρ as the "linear absorption
# coefficient" within these functions (note: throughout the rest of the codebase, we interchangably
# use the simpler term, "absorption coefficient", to refer to this same quantity)

# Combined H I bound-free and free-free opacity
function HI_coefficient(λ, T, Pₑ, H_I_ion_energy=13.598)
    # λ should have units of Ångstroms
    # Pₑ is the partial pressure of the electrons in dyne/cm²
    nₑ = Pₑ / (Korg.kboltz_cgs * T)
    ν = (Korg.c_cgs * 1e8) / λ

    bf_coef = begin
        H_I_partition_val = 2.0 # implicitly in the implementation provided by Gray (2005)
        nH_I = nₑ * 100.0 # this is totally arbitrary
        # pass ne = nHe = 0 to avoid MHD effects. Pass nH = 1 and multiply after the fact for the same reason.
        bf_linear_absorption_coef = Korg.ContinuumAbsorption.H_I_bf([ν], T, 1.0, 0, 0,
                                                                    1 / H_I_partition_val)[1] * nH_I
        bf_linear_absorption_coef / (Pₑ * nH_I)
    end

    ff_coef = begin
        nH_total = nₑ * 100.0 # this is totally arbitrary
        ρ = nH_total * 1.67e-24 / 0.76 # this is totally arbitrary

        #Assume U_I(T) = 2.0 and U_II(T) = 1.0, as in Gray (2005)
        Us = Dict([Korg.species"H_I" => (T -> 2.0), Korg.species"H_II" => (T -> 1.0)])
        χs = [(Korg.RydbergH_eV, -1.0, -1.0)]
        wII, wIII = Korg.saha_ion_weights(T, nₑ, 1, χs, Us)

        nH_I = nH_total / (1 + wII)
        nH_II = nH_total * wII / (1 + wII)

        # compute the linear absorption coefficient  = dτ/ds = opacity*ρ
        ff_linear_absorption_coef = [0.0]
        Korg.ContinuumAbsorption.positive_ion_ff_absorption!(ff_linear_absorption_coef, [ν], T,
                                                             Dict([Korg.species"H II" => nH_II]),
                                                             nₑ)
        ff_linear_absorption_coef[1] / (Pₑ * nH_I)
    end

    bf_coef + ff_coef
end

function Hminus_bf_coefficient(λ, T, Pₑ, ion_energy_H⁻=0.7552)
    # λ should have units of Ångstroms
    # Pₑ is the partial pressure of the electrons in dyne/cm²
    ne = Pₑ / (Korg.kboltz_cgs * T)

    # the values of ndens, ρ and nH_I shouldn't actually matter since we will
    # just divide out all dependence on these variables.
    nH_I, nH, ρ = _semi_realisitic_dens(ne)

    partition_func = 2.0 # may want to include the temperature dependence of the partition function
    ν = (Korg.c_cgs * 1e8) / λ

    linear_absorb_coef = Korg.ContinuumAbsorption.Hminus_bf([ν], T, nH_I / partition_func, ne)[1]
    linear_absorb_coef / (Pₑ * nH_I)
end

function Hminus_ff_coefficient(λ, T, Pₑ)
    # λ should have units of Ångstroms
    # Pₑ is the partial pressure of the electrons in dyne/cm²
    ne = Pₑ / (Korg.kboltz_cgs * T)

    # the values of ndens, ρ and nH_I shouldn't actually matter since we will
    # just divide out all dependence on these variables.
    nH_I, nH, ρ = _semi_realisitic_dens(ne)

    partition_func = 2.0 # may want to include the temperature dependence of the partition function
    ν = (Korg.c_cgs * 1e8) / λ

    linear_absorb_coef = Korg.ContinuumAbsorption.Hminus_ff([ν], T, nH_I / partition_func, ne)[1]
    linear_absorb_coef / (Pₑ * nH_I)
end

# computes the combine H₂⁺ free-free and bound-free absorption in units of cm^2 per H atom (not a
# typo)
# There seems to be an error in our function to compute H₂⁺
function H2plus_coefficient(λ, T, Pₑ)
    # λ should have units of Ångstroms
    # Pₑ is the partial pressure of the electrons in dyne/cm²
    ne = Pₑ / (Korg.kboltz_cgs * T)

    # in the future, we may want to use the abundances, ionization_energies and partition functions
    # used in Gray (2005) since we are comparing to his plots
    ne_div_nH = free_electrons_per_Hydrogen_particle(ne, T)
    nH = ne / ne_div_nH

    wII, wIII = Korg.saha_ion_weights(T, ne, 0x01, Korg.ionization_energies,
                                      Korg.default_partition_funcs)
    nH_I = nH / (1 + wII)
    nH_II = nH * wII / (1 + wII)

    ν = (Korg.c_cgs * 1e8) / λ
    linear_absorb_coef = Korg.ContinuumAbsorption.H2plus_bf_and_ff([ν], T, nH_I, nH_II)[1]
    linear_absorb_coef / (Pₑ * nH_I)
end

# compute He⁻ free-free absorption in units of cm^2 per H atom (not a typo)
function Heminus_ff_coefficient(λ, T, Pₑ)
    # λ should have units of Ångstroms
    # Pₑ is the partial pressure of the electrons in dyne/cm²
    ne = Pₑ / (Korg.kboltz_cgs * T)

    # the values of ndens, ρ and nH_I shouldn't actually matter since we will
    # just divide out all dependence on these variables.

    # In the interest of semi-realistic numbers, Let's assume all free
    # electrons come from Hydrogen
    nH_I, nH, ρ = _semi_realisitic_dens(ne)

    # we need to compute the number density of He I. Since this function exists to facillitate
    # comparisons against plots from Gray (2005), we should employ the abundances from
    # their book. They record that there is 8.51e-2 He particles per H particle in table 16.3
    nHe = nH * 8.51e-2
    # Gray accounts for only neutral and singly ionized helium in his calculation.
    # We use the ionization energy tabulated in table D.1 (χ(He I) = 24.587).
    # Table D.2 of Gray (2005) records the partion function of He II to be
    # as 10^0.301, the closest value at the table's precision to 2.0
    UI, UII = (1.0, 10^-0.301)
    # Because this neglects He III, we don't use Korg.saha_ion_weights.
    wII = (2.0 / ne * (UII / UI) * Korg.translational_U(Korg.electron_mass_cgs, T)
           * exp(-24.587 / (Korg.kboltz_eV * T)))
    nHe_I = nHe * 1 / (1 + wII)

    ν = (Korg.c_cgs * 1e8) / λ
    linear_absorb_coef = Korg.ContinuumAbsorption.Heminus_ff([ν], T, nHe_I / UI, ne)[1]
    linear_absorb_coef / (Pₑ * nH_I)
end

# There appears to be some errors in the H₂⁺ opacities, skipping them for now
const Gray05_opacity_form_funcs = Dict("H" => (HI_coefficient,
                                               Interval(0, Inf), "H I bf and ff"),
                                       "Hminus_bf" => (Hminus_bf_coefficient,
                                                       closed_interval(2250.0, 15000.0),
                                                       "H⁻ bound-free"),
                                       "Hminus_ff" => (Hminus_ff_coefficient,
                                                       closed_interval(2604.0, 113918.0),
                                                       "H⁻ free-free"),
                                       "Heminus_ff" => (Heminus_ff_coefficient,
                                                        closed_interval(5063.0, 151878.0),
                                                        "He⁻ free-free"),
                                       "H2plus" => (H2plus_coefficient,
                                                    closed_interval(3847.0, 25000.0),
                                                    "H₂⁺ ff and bf"))

# this function is defined to make tests easier
function Gray05_comparison_vals(panel, opacity_func_name)
    panel_prop, panel_dict = load_panel_data(panel)
    orig_λ_vals, ref_data = panel_dict[opacity_func_name]
    temperature = panel_prop["temperature"]
    Pₑ = 10^panel_prop["log10(electron_pressure)"]

    func, bounds = Gray05_opacity_form_funcs[opacity_func_name][1:2]

    w = contained.(orig_λ_vals, Ref(bounds))
    calculated_vals = func.(orig_λ_vals[w], temperature, Pₑ)
    (calculated_vals * 1e26, ref_data[w])
end

end # end the definition of the Gray_opac_compare module

# ================================================================================================

# Below, we define a set of functions that are used to compare our implementation of bound-free
# opacity calculations using data from the Opacity Project against calculations of the same
# opacities using the hydrogenic aproximation. This nominally means that we're comparing
# opacities from H I and He II
#
# In the future, it might also be nice to include comparisons for multi-electron data against
# experimental data.
module OP_compare
using Korg

function calc_hydrogenic_bf_absorption_coef(λ_vals, T, ndens_species, spec; use_OP_data=true)
    @assert spec in Korg.Species.(["H I", "He II"])
    @assert (spec == Korg.species"H I") || use_OP_data

    νs = Korg.c_cgs ./ (λ_vals * 1e-8)
    if use_OP_data || spec == Korg.species"He II"
        σ_itp = Korg.ContinuumAbsorption.metal_bf_cross_sections[spec]
        exp.(log(ndens_species) .+ σ_itp.(νs, log10(T))) * 1e-18 #convert to cm^2
    else # spec is H I
        invU = 1 / Korg.default_partition_funcs[spec](log(T))
        # assume nHe = ne = 0 for the purposed of MHD
        reverse(Korg.ContinuumAbsorption.H_I_bf(reverse(νs), T, ndens_species, 0, 0, invU))
    end
end

_dflt_λ_vals(spec) =
    if spec == Korg.species"H I"
        # skip the breaks because that's where MHD make a difference.  Don't go to far red because
        # those levels get dissolved.
        [500:1.0:600; 3000:1.0:3100; 4000:1.0:4100; 10_000:1.0:10_100]
    else # He II
        500:1.0:900 # cross-section is below Float32 min above ~912 Å
    end

_hydrogenic_rtol(species_name) = (species_name == "H_I") ? 0.12 : 0.30
end # end the definition of OP_compare module
