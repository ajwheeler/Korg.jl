using Base
using Interpolations

# While the StateID struct and its helper functions may seem a little over-engineered, they
# primarily exist to document the iSLP and iLV quantities used by TOPBase

"""
    StateID(iSLP, iLV)

`StateID` represents a given state (electron configuration) of an ion species. The arguments are
taken from the TOPBase.

# Arguments
- `iSLP::UInt16`: An integer representing the quantum atomic numbers SLπ (total spin, orbital
  angular momentum, and parity). The value is given by `100*(2*S+1) + 10*L + π`. Note that TOPBase
  documentation refers to the parity quantum number as P or π, but for consistency with how TOPBase
  presents the data, we prefer P in our implementation code.
- `iLV::UInt8`: Distinguishes between configurations that have the same values of `iSLP`. TOPBase
  seems to index the configuration based on the energy level of the configuration (with `iLV = 1`
  corresponding to the lowest energy).

# Notes
The combination of this struct and an ion species (or neutral element) should be enough to fully
specify a given electron configuration. For example, consider He I:
- `iSLP = 100, iLV = 1` has the electron configuration: 1s²
- `iSLP = 100, iLV = 2` has the electron configuration: 1s 2s
- `iSLP = 100, iLV = 3` has the electron configuration: 1s 3s
- `iSLP = 111, iLV = 1` has the electron configuration: 1s 2p
- `iSLP = 111, iLV = 2` has the electron configuration: 1s 3p
"""
struct StateID
    iSLP::UInt16
    iLV::UInt8
    function StateID(iSLP, iLV)
        @assert 99 < iSLP < 992 # this upper limit is set based on the data source
        @assert (iSLP % 10) < 2 # trailing digit must be 0 or 1
        @assert 0 < iLV # the upper limit is set based on the data source
        new(iSLP, iLV)
    end
end

Base.:(<)(a::StateID, b::StateID) = (a.iSLP < b.iSLP) || ((a.iSLP == b.iSLP) && (a.iLV < b.iLV))
Base.show(io::IO, obj::StateID) = print(io, "StateID{iSLP = ", obj.iSLP, ", iLV = ", obj.iLV, "}")

_parity_vals = ("even", "odd") # predefine these to avoid allocations

function _quantum_prop_type(state_id::StateID)
    S_code = div(state_id.iSLP, 100) # hundreds place of state_id.iSLP
    twice_S = S_code - 1
    L = div(state_id.iSLP - S_code*100, 10) # tens place of state_id.iSLP
    parity = _parity_vals[(state_id.iSLP % 2) + 1] # 0 or 1 in ones place of state_id.iSLP
    (Int32(twice_S), Int32(L), parity)
end
get_S_quantum_num(state_id::StateID) = _quantum_prop_type(state_id)[1] / 2.0f0
get_L_quantum_num(state_id::StateID) = _quantum_prop_type(state_id)[2]
get_parity(state_id::StateID) = _quantum_prop_type(state_id)[3]

function get_statistical_weight(state_id::StateID)
    twice_S, L, _ = _quantum_prop_type(state_id)
    (twice_S+1)*(2*L+1)
end

"""
    StateProp(ion_energy_ryd, excitation_potential_ryd, statistical_weight, electron_config)

`StateProp` holds data associated with different electron configurations. It's used to hold
data read from the TOPBase cross-section queries and to compute the statistical weight.
"""
struct StateProp
    ion_energy_ryd::Float64 # energy with respect to ionization potential in Rydbergs
    excitation_potential_ryd::Float64 # energy above ground state in Rydbergs
    statistical_weight::UInt16
    electron_config::String
end


# define some generally useful functions related to parsing
function _parse_species_stateid(str::AbstractString, start_Z::Integer,
                                start_numElectrons::Integer, start_iSLP::Integer,
                                start_iLV::Integer, last_iLV::Integer)
    # parses the atomic_num, number of electrons, iSLP, and iLV entries of the table
    atomic_num = parse(UInt8, str[start_Z:start_numElectrons-1])
    num_electrons = parse(UInt8, str[start_numElectrons:start_iSLP-1])
    iSLP = parse(UInt16, str[start_iSLP:start_iLV-1])
    iLV = parse(UInt16, str[start_iLV:last_iLV])
    (Species(Formula(atomic_num), atomic_num - num_electrons), StateID(iSLP, iLV))
end

# Functionallity for reading photoionization cross sections. Since these tables can be massive, we
# use an iterator to do this. This is based on the design used in the Julia standard library for
# implementing the ``EachLine`` iterable object. That code used the following license:
#
# Copyright (c) 2009-2021: Jeff Bezanson, Stefan Karpinski, Viral B. Shah,
# and other contributors:
#
# https://github.com/JuliaLang/julia/contributors
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


const _TABLE_BREAK = " ================================================"
const _HEADER_LINE_LENGTH = length(_TABLE_BREAK)

struct EachPhotoIonSubtable{IOT <: IO}
    stream::IOT
    ondone::Function
end

"""
    each_photo_ion_subtable(filename::AbstractString)

Create an iterable `EachPhotoIonSubtable` object for parsing photoionization cross-sections.

The input file is organized as a series of subtables. Each subtable provides the photoionization
energies (make sure this isn't the energy of the photoionized electron) in Ryd, and the associated
cross section (in megabarns) for a different state.

Each entry in the `EachPhotoIonSubtable` iterable provides a named tuple with the keys:
- name: holds the species name
- state_id: holds a `StateID` instance
- photon_energy_ryd: the photoionization values (check), in Ryd, stored in an array of Float32 vals
- cross_section_MegaBarn: the associated photoionization cross-sections. This array should have the
  same shape and datatype as photon_energy_ryd

The underlying file is only closed after the `EachPhotoIonSubtable` is garbage collected.

#Notes
This implementation is not robust enough to deal with blank lines at the start or end of the file.

The table should have been produced by a
[TOPbase Photoionisation Cross Section query](http://cdsweb.u-strasbg.fr/topbase/xsections.html)
that uses the following query options:
- Only a single atomic number and a single number of electrons is included in the query (other
  choices may yield too many results, which causes TOPBase to truncate output).
- There isn't any limitation on the number the quantum numbers, the energy or the levels
- The query request levels ordered in the "Level order in each series"
"""
function each_photo_ion_subtable(filename::AbstractString)
    stream = open(filename, "r")
    # confirm that standard file header is present and move the position of the
    # stream past this point
    @assert readline(stream) == _TABLE_BREAK
    @assert readline(stream) == "       I  NZ  NE  ISLP  ILV        E(RYD)      NP"
    @assert readline(stream) == _TABLE_BREAK
    EachPhotoIonSubtable(stream, ()->close(stream))::EachPhotoIonSubtable
end

function Base.iterate(itr::EachPhotoIonSubtable, state=nothing)
    next_line = readline(itr.stream)

    if next_line == _TABLE_BREAK
        return (itr.ondone(); nothing)
    end

    # first parse the subtable's header
    species, state_id = _parse_species_stateid(next_line, 9, 13, 17, 23, 27)
    ion_energy_ryd = parse(Float32, next_line[28:41])

    num_entries = parse(UInt32, next_line[42:49])

    # this could probably be 32-bit.
    photon_energy_ryd = Array{Float32,1}(undef, num_entries)
    cross_section_MegaBarn = Array{Float32,1}(undef, num_entries)
    for i = 1:num_entries # this properly handles the case when num_entries == 0
        cur_line = readline(itr.stream)
        photon_energy_ryd[i] = parse(Float32, cur_line[1:14])
        cross_section_MegaBarn[i] = parse(Float32, cur_line[15:24])
    end

    # parse the species, and the level
    ((species = species, state_id=state_id, ion_energy_ryd = ion_energy_ryd,
      photon_energy_ryd_arr = photon_energy_ryd,
      cross_section_MegaBarn_arr = cross_section_MegaBarn), nothing)
end

# this and the following definition implement Julia's informal iterator interface
Base.eltype(::Type{<:EachPhotoIonSubtable}) =
    Tuple{String, StateID, Float32, Array{Float32,1}, Array{Float32,1}}
Base.IteratorSize(::Type{<:EachPhotoIonSubtable}) = SizeUnknown()

# In a couple tables, there's a place where the photon energy is listed out of order. This seems to
# be okay given that this is in a section of the table where the cross-section is zero.
const _species_with_problematic_subtables = (Species("He_I"), Species("C_I"), Species("C_II"),
                                             Species("Al_II"), Species("O_II"))

function _interpolator_from_subtable(table_photon_energy_ryd::Vector{Float32},
                                     table_cross_section_MBarn::Vector{Float32},
                                     species::Union{Species,Nothing} = nothing,
                                     extrapolation_bc = 0.0)
    if species in _species_with_problematic_subtables
        
        for j in 1:length(table_photon_energy_ryd)-1
            if table_photon_energy_ryd[j+1] < table_photon_energy_ryd[j]
                @assert j > 1
                @assert (j + 2) <= length(table_photon_energy_ryd)
                @assert table_photon_energy_ryd[j+2] >= table_photon_energy_ryd[j+1]
                @assert table_photon_energy_ryd[j+2] >= table_photon_energy_ryd[j]
                @assert table_photon_energy_ryd[j+1] >= table_photon_energy_ryd[j-1]
                @assert table_photon_energy_ryd[j] >= table_photon_energy_ryd[j-1]

                @assert table_cross_section_MBarn[j+2] == 0.0
                @assert table_cross_section_MBarn[j+1] == 0.0
                @assert table_cross_section_MBarn[j] == 0.0
                @assert table_cross_section_MBarn[j-1] == 0.0

                # swap the order of the 2 values
                temp = table_photon_energy_ryd[j]
                table_photon_energy_ryd[j] = table_photon_energy_ryd[j+1]
                table_photon_energy_ryd[j+1] = temp
            end
        end
    end
    LinearInterpolation(table_photon_energy_ryd, table_cross_section_MBarn,
                        extrapolation_bc=extrapolation_bc)
end

function _try_env_prefix_path(prefix_env_var, fname)
    prefix = get(ENV, prefix_env_var, nothing)
    if isnothing(prefix)
        e = ErrorException(string("The \"", prefix_env_var, "\" environment variable was not set"))
        throw(e);
    end
    path = joinpath(prefix, fname)
    if !isfile(path) # this properly handles the case where path is a symlink to a file
        throw(ErrorException(string("There is no file at ", path)));
    end
    path
end

# return the minimum iterator over cross_section subtables and the ionization energy of the ground
# state (according to the table)
function _get_cross_sec_data(species::Species,
                             cross_sec_file::Union{AbstractString,Nothing} = nothing)

    # get the photo-ionization subtable iterator:
    my_cross_sec_file = if isnothing(cross_sec_file)
        # determine the string representation of Species
        # we don't directly use repr(species) since the representation could change in the future
        # to not use roman numeral (e.g. if support is added for negatively charged species)
        formula_str = repr("text/plain", species.formula)
        roman_numeral = get_roman_numeral(species)
        fname = string(formula_str, "_", roman_numeral, ".txt")

        println("cross_sec_file isn't specified. Searching for file in the directory specified ",
                "by the KORG_OP_DIR environment variable.")
        _try_env_prefix_path("KORG_OP_DIR", fname)
    else
        cross_sec_file
    end

    # determine the ionization energy of the ground state (from the subtables)
    gs_ion_energy_ryd = begin # this could be optimized significantly!
        cur_min = floatmax(Float32)
        for (cur_species, _, ion_energy_ryd, _, _) in each_photo_ion_subtable(my_cross_sec_file)
            @assert (cur_species == species) # this will catch mistakes in the query that generated
                                             # the cross_section data file 
            cur_min = min(ion_energy_ryd, cur_min)
        end
        cur_min
    end

    subtable_iter = each_photo_ion_subtable(my_cross_sec_file)
    subtable_iter, gs_ion_energy_ryd
end

@doc raw"""
    weighted_bf_cross_section_TOPBase(λ_vals, T_vals, species; extrapolation_bc = 0.0
                                      cross_sec_file = nothing,
                                      partition_func = partition_funcs[species])

Calculate the combined bound-free cross section (in cm²) of an ion species using data from the
Opacity Project. The result is the mean of the the bound-free cross section of each electron state,
weighted by their Boltzmann probabilities. It includes the LTE correction for stimulated emission. 
The result is stored in a 2D array with a size given by `(length(λ_vals), length(T_vals))`.

# Arguments
- `λ_vals::AbstractVector`: wavelengths (in cm) to compute the cross sections at
- `T_vals::AbstractVector`: temperatures (in K) to compute the cross sections at
- `species::Species`: the starting atom or atomic ion of the bound-free reaction.

# Keyword Arguments
- `extrapolation_bc`: Specifies handling of extrapolation during linear interpolation (this is
  passed to `Interpolations.LinearInterpolation`).
- `cross_sec_file`: Optional keyword that specifies the path to the text file that holds the
  cross-section subtables. When `nothing` (the default), this function tries to load the file
  from `$KORG_OP_DIR/<atomic-symbol>_<ion-state>.txt`. `<atomic-symbol>` is replaced with the
  atomic symbol of the species (e.g. `H`, `He`, `Al`) and `<ion-state>` is replaced with the 
  capital Roman numeral(s) specifying the ionization state.
- `partition_func`: Specifies the partition function for the current species. This should be a
  callable that accepts temperature as an argument. When this is `nothing`, it falls back to the
  default partition functions.

# Explanation
In more mathematical rigor, this function basically evaluates the following equation:

``\overline{\sigma_{\rm bf}}(\lambda,T) = \sum_{{\bf n}\in{\rm states}} \frac{ g_{\bf n}}{U(T)} e^{-E_{\bf n}/(k_B T)}\left(1 - e^{-h c/(\lambda k_B T)}\right) \sigma_{{\rm bf}, {\bf n}}(\lambda)``

In this equation:
- ``g_{\bf n}`` is the statistical weight of state ``{\bf n}``.
- ``U`` is the partition function.
- ``E_{\bf n}`` is the excitation potential of state ``{\bf n}``. In other words, its the energy 
  difference between state ``{\bf n}`` and the ground state.
- ``\sigma_{\lambda,{\rm bf}, {\bf n}}`` is the bound-free cross-section (a.k.a. photo-ionization
  cross section) of state ``{\bf n}``.

The linear aborption coefficient for bound-free absorption is simply:
``\alpha_{\lambda,{\rm bf}} = n_{\rm tot} \overline{\sigma_{\rm bf}}(\lambda,T)``.

# Notes
In the future, it needs to be possible to pass an argument that adjusts the energy levels of
different electron configurations.
"""
function weighted_bf_cross_section_TOPBase(λ_vals::AbstractVector, T_vals::AbstractVector,
                                           species::Species; extrapolation_bc=0.0,
                                           cross_sec_file::Union{AbstractString,Nothing} = nothing,
                                           partition_func = partition_funcs[species])
    # Note: it would not take a lot of work to make this work with 32-bit floats (since the tables
    # only store 32-bit floats). If doing that, it might be wise to accept λ in units of ångstroms
    # (for precision purposes)

    # TODO: make it possible to pass in a dict mapping dict names to the various energy levels
    #       (to allow for overwriting the values in the table with empirical values)

    @assert !ismolecule(species) # species must be an atom or atomic ion

    # precompute Temperature-dependent constant
    inv_partition_func_val = 1.0./partition_func.(T_vals)

    β_Ryd = RydbergH_eV./(kboltz_eV .* T_vals)

    # convert λ_vals to photon energies
    photon_energies = (hplanck_eV * c_cgs / RydbergH_eV) ./ λ_vals

    # prepare the output array where results will be accumulated
    weighted_average = zeros(eltype(photon_energies), (length(photon_energies),length(T_vals)))

    cross_section_subtables, gs_ion_energy_ryd = _get_cross_sec_data(species, cross_sec_file)

    for subtable in cross_section_subtables
        if length(subtable.photon_energy_ryd_arr) == 0
            continue
        end

        cross_section_interp = _interpolator_from_subtable(subtable.photon_energy_ryd_arr,
                                                           subtable.cross_section_MegaBarn_arr,
                                                           species, extrapolation_bc)
        statistical_weight = get_statistical_weight(subtable.state_id)
        excitation_potential_ryd = subtable.ion_energy_ryd - gs_ion_energy_ryd

        # now iterate over Temperatures
        for j in 1:length(T_vals)

            # now compute the weighting of the current state under LTE
            weight = inv_partition_func_val[j] * statistical_weight * exp(-excitation_potential_ryd
                                                                          * β_Ryd[j])

            current_cross_section = ( cross_section_interp.(photon_energies) .*
                                      (1.0 .- exp.(-photon_energies.*β_Ryd[j])) )
            view(weighted_average, :, j) .+= (current_cross_section .* weight)
        end
    end

    # convert weighted_average units from megabarn to cm²
    weighted_average .*= 1e-18
    weighted_average
end


"""
    absorption_coef_bf_TOPBase(λ, T, ndens_species, species; kwargs...)

Computes the bound-free linear absorption coefficient, α, (in units of cm⁻¹) using data from the 
Opacity Project (downloaded from TOPBase) and assuming LTE. The result is stored in a 2D array with
a size given by `(length(λ), length(T))`.

# Arguments
- `λ::AbstractVector`: wavelengths (in cm) to compute α at.
- `T::AbstractVector`: temperatures (in K) to compute α at.
- `ndens_species::AbstractVector`: number density (in cm⁻³) of the species at each temperature 
  (this should have the same length as T).
- `species::Species`: the starting atom or atomic ion of the bound-free reaction.

# Keyword Arguments
This function accepts all the same keyword arguments as `weighted_bf_cross_section_TOPBase`.
"""
function absorption_coef_bf_TOPBase(λ::AbstractVector, T::AbstractVector,
                                    ndens_species::AbstractVector, species::Species; kwargs...)
    @assert size(T) == size(ndens_species)

    arr = weighted_bf_cross_section_TOPBase(λ, T, species; kwargs...)

    # to save time and memory, let's modify the values of arr in-place
    for (i,n) in enumerate(ndens_species)
        view(arr,:,i) .*= n
    end

    arr
end
