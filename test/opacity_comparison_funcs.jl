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

# compute the ratio of the nₑ to n_H (the number density of all H I and H II)
# this function could be reused in the future for testing solutions to coupled Saha equations
function free_electrons_per_Hydrogen_particle(nₑ, T, abundances = Korg.solar_abundances,
                                              ionization_energies = Korg.ionization_energies,
                                              partition_funcs = Korg.partition_funcs)
    out = 0.0
    for element in 1:Korg.Natoms
        wII, wIII = Korg.saha_ion_weights(T, nₑ, element, Korg.ionization_energies, 
                                          Korg.partition_funcs)

        nₑ_per_ndens_species = wII/(1 + wII + wIII) + 2wIII/(1 + wII + wIII)
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
function _semi_realisitic_dens(ne::F, fion::F= 0.02, HydrogenMassFrac::F= 0.76) where F <: Real
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
    nₑ =  Pₑ/(Korg.kboltz_cgs * T)
    ν = (Korg.c_cgs*1e8)/λ

    bf_coef = begin
        H_I_partition_val = 2.0 # implicitly in the implementation provided by Gray (2005)
        nH_I = nₑ * 100.0 # this is totally arbitrary
        ρ = nH_I * 1.67e-24/0.76 # this is totally arbitrary
        bf_opac = Korg.ContinuumOpacity.H_I_bf(nH_I/H_I_partition_val, ν, ρ, T, H_I_ion_energy)
        bf_opac * ρ / (Pₑ * nH_I)
    end

    ff_coef = begin
        nH_total = nₑ * 100.0 # this is totally arbitrary
        ρ = nH_total * 1.67e-24/0.76 # this is totally arbitrary
        
        #Assume U_I(T) = 2.0 and U_II(T) = 1.0, as in Gray (2005)
        Us = Dict([Korg.literals.H_I=>(T -> 2.0), Korg.literals.H_II=>(T -> 1.0)])
        χs = [(Korg.RydbergH_eV, -1.0, -1.0)]
        wII, wIII = Korg.saha_ion_weights(T, nₑ, 1, χs, Us)

        nH_I = nH_total / (1 + wII)
        nH_II = nH_total * wII/(1 + wII)
        ff_opac = Korg.ContinuumOpacity.H_I_ff(nH_II, nₑ, ν, ρ, T)
        ff_opac * ρ / (Pₑ * nH_I)
    end

    bf_coef + ff_coef
end

function Hminus_bf_coefficient(λ, T, Pₑ, ion_energy_H⁻ = 0.7552)
    # λ should have units of Ångstroms
    # Pₑ is the partial pressure of the electrons in dyne/cm²
    ne =  Pₑ/(Korg.kboltz_cgs * T)

    # the values of ndens, ρ and nH_I shouldn't actually matter since we will
    # just divide out all dependence on these variables.
    nH_I, nH, ρ = _semi_realisitic_dens(ne)

    partition_func = 2.0 # may want to include the temperature dependence of the partition function
    ν = (Korg.c_cgs*1e8)/λ

    opacity = Korg.ContinuumOpacity.Hminus_bf(nH_I/partition_func, ne, ν, ρ, T, ion_energy_H⁻)
    opacity * ρ / (Pₑ * nH_I)
end

function Hminus_ff_coefficient(λ, T, Pₑ)
    # λ should have units of Ångstroms
    # Pₑ is the partial pressure of the electrons in dyne/cm²
    ne =  Pₑ/(Korg.kboltz_cgs * T)

    # the values of ndens, ρ and nH_I shouldn't actually matter since we will
    # just divide out all dependence on these variables.
    nH_I, nH, ρ = _semi_realisitic_dens(ne)

    partition_func = 2.0 # may want to include the temperature dependence of the partition function
    ν = (Korg.c_cgs*1e8)/λ

    opacity = Korg.ContinuumOpacity.Hminus_ff(nH_I/partition_func, ne, ν, ρ, T)
    opacity * ρ / (Pₑ * nH_I)
end

# computes the combine H₂⁺ free-free and bound-free absorption in units of cm^2 per H atom (not a
# typo)
# There seems to be an error in our function to compute H₂⁺
function H2plus_coefficient(λ, T, Pₑ)
    # λ should have units of Ångstroms
    # Pₑ is the partial pressure of the electrons in dyne/cm²
    ne =  Pₑ/(Korg.kboltz_cgs * T)

    # in the future, we may want to use the abundances, ionization_energies and partition functions
    # used in Gray (2005) since we are comparing to his plots
    ne_div_nH = 
        free_electrons_per_Hydrogen_particle(ne, T)
    nH = ne/ne_div_nH

    wII, wIII = Korg.saha_ion_weights(T, ne, 0x01, Korg.ionization_energies, 
                                         Korg.partition_funcs)
    nH_I = nH / (1 + wII)
    nH_II = nH * wII / (1 + wII)


    # set the partition function to 2.0 in nH_I_div_partition so that n(H I, n=1) = n(H I) for this
    # calculation, which is an assumption that Gray (2005) implicitly uses
    nH_I_div_partition = nH_I/2.0
    ρ = 1.0 # arbitrary value because we divide it out after

    ν = (Korg.c_cgs*1e8)/λ
    opacity = Korg.ContinuumOpacity.H2plus_bf_and_ff(nH_I_div_partition, nH_II, ν, ρ, T)
    opacity * ρ / (Pₑ * nH_I)
end



# compute He⁻ free-free absorption in units of cm^2 per H atom (not a typo)
function Heminus_ff_coefficient(λ, T, Pₑ)
    # λ should have units of Ångstroms
    # Pₑ is the partial pressure of the electrons in dyne/cm²
    ne =  Pₑ/(Korg.kboltz_cgs * T)

    # the values of ndens, ρ and nH_I shouldn't actually matter since we will
    # just divide out all dependence on these variables.

    # In the interest of semi-realistic numbers, Let's assume all free
    # electrons come from Hydrogen
    nH_I, nH, ρ = _semi_realisitic_dens(ne)

    # we need to compute the number density of He I. Since this function exists to facillitate
    # comparisons against plots from Gray (2005), we should employ the abundances from
    # their book. They record that there is 8.51e-2 He particles per H particle in table 16.3
    nHe = nH*8.51e-2
    # Gray accounts for only neutral and singly ionized helium in his calculation.
    # We use the ionization energy tabulated in table D.1 (χ(He I) = 24.587). 
    # Table D.2 of Gray (2005) records the partion function of He II to be
    # as 10^0.301, the closest value at the table's precision to 2.0
    UI, UII = (1.0, 10^-0.301)
    # Because this neglects He III, we don't use Korg.saha_ion_weights.
    wII = (2.0/ne * (UII/UI) * Korg.translational_U(Korg.electron_mass_cgs, T) 
           * exp(-24.587/(Korg.kboltz_eV*T)))
    nHe_I = nHe * 1/(1+wII)
                         
    ν = (Korg.c_cgs*1e8)/λ
    opacity = Korg.ContinuumOpacity.Heminus_ff(nHe_I/UI, ne, ν, ρ, T)
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

function inbounds(bounds::Bounds, vals::Array{F}) where F <: Real
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
# unlike for H₂⁺, this could simply be a consequence of the different approximations that are used
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

function calc_hydrogenic_bf_absorption_coef(λ_vals,  T, ndens_species, species_name;
                                            use_OP_data = false)
    @assert species_name in ["H_I", "He_II"]

    if use_OP_data
        cross_sec_file = if species_name == "H_I"
            joinpath(@__DIR__, "data/TOPbase_cross_section_H_I.txt")
        else
            joinpath(@__DIR__, "data/TOPbase_cross_section_He_II.txt")
        end
        Korg.ContinuumOpacity.absorption_coef_bf_TOPBase(λ_vals.*1e-8, [T], [ndens_species],
                                                         Korg.Species(species_name);
                                                         extrapolation_bc = 0.0,
                                                         cross_sec_file = cross_sec_file)[:, 1]
    else
        ν_vals = (Korg.c_cgs*1e8)./λ_vals # Hz
        ndens_div_partition = ndens_species/Korg.partition_funcs[Korg.Species(species_name)](T)
        ρ = 1.0 # this is unphysical, but it works out fine for a hydrogenic atom
        κ_dflt_approach = if species_name == "H_I"
            H_I_ion_energy = 13.598
            Korg.ContinuumOpacity.H_I_bf.(ndens_div_partition, ν_vals, ρ, T, H_I_ion_energy)
        else
            He_II_ion_energy = 54.418
            Korg.ContinuumOpacity.He_II_bf.(ndens_div_partition, ν_vals, ρ, T, He_II_ion_energy)
        end
        κ_dflt_approach/ρ
    end
end

# the following functions could all be refactored so that they are more concise... (they could
# honestly just be constants)

function _dflt_λ_vals(species_name)
    first, last = (species_name == "H_I") ? (80,80000) : (25, 40000)
    map((λ) -> Float32(λ), first:1.0:last)
end

function _λ_comp_intervals(species_name)
    # because the energy state data used by the Opacity Project is somewhat inaccurate (note the
    # disagreement gets worse at higher energy levels), we can only make meaningful comparisons
    # over specific wavelength intervals.

    if species_name == "H_I"
        [80.0       911.0;
         1000.0     3645.0;
         4000.0     8203.0;
         9300.0     14579.0;
         18000.0    22789.0]
    else
        [25.0       10.0^2.35;
         10.0^2.37  10.0^2.95;
         10.0^2.99  10.0^3.3;
         10.0^3.37  10.0^3.55;
         10.0^3.65  10.0^3.75]
    end
end

_hydrogenic_rtol(species_name) = (species_name == "H_I") ? 0.11 : 0.06
end # end the definition of OP_compare module
