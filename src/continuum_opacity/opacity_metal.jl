# the only reason that we NEED to load in the electron_configurations data is that it's not
# completely clear to me exactly how I would compute certain quantities (i.e. the statistical
# weight and the ground state energy)

using Base
using Interpolations
# import Korg.partition_funcs, Korg.atomic_symbols, Korg.numerals, and constants


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
Base.:(==)(a::StateID, b::StateID) = (a.iSLP == b.iSLP) && (a.iLV == b.iLV)
Base.hash(obj::StateID) = Base.hash((UInt32(obj.iSLP) << 16) + obj.iLV)
Base.show(io::IO, obj::StateID) = print(io, "StateID{iSLP = ", obj.iSLP, ", iLV = ", obj.iLV, "}")

function _quantum_prop_type(state_id::StateID)
    str = String(iSLP)
    S = div(Integer(str[1])-1,2)
    L = Integer(str[2])
    parity = ["even", "odd"][Integer(str[3])+1]
    (S, L, parity)
end
get_S_quantum_num(state_id::StateID) = _quantum_prop_type(state_id)[1]
get_L_quantum_num(state_id::StateID) = _quantum_prop_type(state_id)[2]
get_parity(state_id::StateID) = _quantum_prop_type(state_id)[3]


"""
    StateProp(ion_energy_ryd, excitation_potential_ryd, statistical_weight, electron_config)

`StateProp` holds data associated with different electron configurations. It's mostly used to hold
data read from TOPBase electron configuration queries (see `read_electron_configurations`)

# Notes

It should be possible to directly determine `electron_config` and `statistical_weight` given a
species name and a `StateID` instance. However, the procedure for doing so is not entirely obvious
"""
struct StateProp
    ion_energy_ryd::Float64 # energy with respect to ionization potential in Rydbergs
    excitation_potential_ryd::Float64 # energy above ground state in Rydbergs
    statistical_weight::UInt16
    electron_config::String
end


# define some generally useful functions related to parsing
function _parse_Z_numelectrons_stateid(str::AbstractString, start_Z::Integer,
                                       start_numElectrons::Integer, start_iSLP::Integer,
                                       start_iLV::Integer, last_iLV::Integer)
    # parses the atomic_num, number of electrons, iSLP, and iLV entries of the table
    atomic_num = parse(UInt8, str[start_Z:start_numElectrons-1])
    num_electrons = parse(UInt8, str[start_numElectrons:start_iSLP-1])
    iSLP = parse(UInt16, str[start_iSLP:start_iLV-1])
    iLV = parse(UInt16, str[start_iLV:last_iLV])
    (atomic_num, num_electrons, StateID(iSLP, iLV))
end

_get_species_string(atomic_num::Integer, num_electrons::Integer) =
    atomic_symbols[atomic_num] * "_" * roman_numerals[atomic_num-num_electrons+1]


# Functionallity for reading tables of electron configuration data
const _ELECTRON_CONFIG_TABLE_BREAK =
        " ===================================================================="

"""
    read_electron_configurations(fname)

Parse the electron configuration table stored in fname. This function returns a vector of 2-tuples.
The first entry of each tuple is the name of a species and the second entry is a dictionary holding
properties for various electron configurations of that species. The keys of each dictionary are
instances of `StateID` and the values are instances of `StateProp`.

# Notes
This implementation is not robust enough to deal with blank lines at the start or end of the file

The table should have been produced by a
[TOPbase Energy levels query](http://cdsweb.u-strasbg.fr/topbase/energy.html). The
requirements for the query are described below:

This function expects a table produced with the following query options:
- Only a single atomic number is included in the query
- There isn't any limitation on the number of electrons, the quantum numbers, the energy or the
  levels
- The query request levels ordered in the "Level order in each series"

Under the output options, the following output data options should be checked:
- "Electron configuration"
- "Energy (Ryd) wrt ionization potential"
- "Energy (Ryd) wrt ground state"
- "Statistical weight"

A lot of the data read in by this function is technically unnecessary:
- The "Electron Configuration" data is only used for debugging purposes.
- It should be possible to construct the "Statistical Weight" data (and "Statistical Weight" data)
  given the species names and a `StateID`. However, the procedure for doing so is not obvious.
- Given "Energy (Ryd) wrt ionization potential" for all states of a given species, we should be
  able to compute "Energy (Ryd) wrt ground state". However, it's not completely obvious how to
  identify the ground state.
"""
function read_electron_configurations(fname)

    # The choice to use a dict for pair StateID instances with StateProp instances was made for
    # simplicity. A more efficient data structure could be chosen since the values are ordered (and
    # are generally accessed in order)

    outputs = Vector{Tuple{String,Dict{StateID,StateProp}}}()

    linelist = open(fname) do f
        file_atomic_num = nothing
        current_num_electrons = nothing
        current_dict = nothing

        # confirm that standard file header is present and move the position of the stream past
        # this point
        @assert readline(f) == _ELECTRON_CONFIG_TABLE_BREAK
        @assert readline(f) ==
            "      i NZ NE iSLP iLV iCONF                 E(RYD)      TE(RYD)   gi"
        @assert readline(f) == _ELECTRON_CONFIG_TABLE_BREAK

        for line in eachline(f)
            if line == _ELECTRON_CONFIG_TABLE_BREAK # the last line
                break
            end

            atomic_num, num_electrons, state_id =
                _parse_Z_numelectrons_stateid(line, 9, 11, 14, 19, 22)
            species = _get_species_string(atomic_num, num_electrons)

            if file_atomic_num === nothing
                file_atomic_num = atomic_num
            else
                @assert file_atomic_num == atomic_num
            end

            if current_num_electrons != num_electrons
                current_dict = Dict{StateID,StateProp}()
                push!(outputs, (species, current_dict))
                current_num_electrons = num_electrons
            end

            electron_config = strip(line[24:39])
            ion_energy_ryd = parse(Float32, line[40:51])
            excitation_potential_ryd = parse(Float32, line[53:64])

            # finally, parse statistical weight. This is always listed as 1 or 2 digits, a decimal
            # point, and then a zero. Since statistical weight is always an int, we'll only parse
            # the leading 2 characters.

            if line[66:69] == "****"
                # I think that this is only an issue for 1 highly excited state of Fe
                continue
            end
            statistical_weight = parse(UInt16, line[66:67])

            current_dict[state_id] = StateProp(ion_energy_ryd, excitation_potential_ryd,
                                               statistical_weight, electron_config)
        end
    end
    outputs
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
    atomic_num, num_electrons, state_id =
        _parse_Z_numelectrons_stateid(next_line, 9, 13, 17, 23, 27)
    species = _get_species_string(atomic_num, num_electrons)
    # ion_energy_ryd is duplicated from the other table. Technically, we don't need the other
    # table. we just need to implement a function to compute the statistical weight from state_id
    # and the species
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
    ((species = species, state_id=state_id, photon_energy_ryd = photon_energy_ryd,
      cross_section_MegaBarn = cross_section_MegaBarn), nothing)
end

Base.eltype(::Type{<:EachPhotoIonSubtable}) =
    Tuple{String, StateID, Array{Float32,1}, Array{Float32,1}}

Base.IteratorSize(::Type{<:EachPhotoIonSubtable}) = SizeUnknown()



function _interpolator_from_subtable(table_photon_energy_ryd::Vector{Float32},
                                     table_cross_section_MBarn::Vector{Float32},
                                     species_name::Union{AbstractString,Nothing} = nothing,
                                     extrapolation_bc = 0.0)
    if species_name in ["He_I", "C_I", "C_II", "Al_II", "O_II"]
        # In a couple tables, there's a place where the photon energy is listed out of order
        # This seems to be okay given that this is in a section of the table where the
        # cross-section is zero.
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

function _get_tabulated_data(species_name, elec_conf_table = nothing,
                             cross_sec_file::Union{AbstractString,Nothing} = nothing)

    # Split appart the species_name (and check that its sensible)
    split_name = split(species_name, "_")
    @assert length(split_name) == 2
    element_name = split_name[1]
    ion_state = split_name[2]
    @assert ion_state in roman_numerals

    # find electron configuration dict:
    my_elec_config_file = if isnothing(elec_conf_table)
        # look at environment variable KORG_OP_ELECTRON_CONFIG_DIR
        println("Searching for electron_config_file in KORG_OP_ELECTRON_CONFIG_DIR")
        joinpath(_data_dir,
                 string("TOPbase/electron_config/", element_name, ".txt"))
    else
        elec_conf_table
    end
    println(my_elec_config_file)

    result = read_electron_configurations(my_elec_config_file)
    rslt_index = 0
    for i = 1:length(result)
        if result[i][1] == species_name
            rslt_index = i
        end
    end
    @assert rslt_index != 0
    data_dict = result[rslt_index][2]

    # get the photo-ionization subtable iterator:
    my_cross_sec_file = if isnothing(cross_sec_file)
        joinpath(_data_dir,
                 string("TOPbase/cross_sections/", species_name, ".txt"))
    else
        cross_sec_file
    end

    (data_dict, each_photo_ion_subtable(my_cross_sec_file))
end


@doc raw"""
    weighted_bf_cross_section_TOPBase(λ_vals, T_vals, species_name; convert_to_cm2 = false,
                                      partition_func = nothing, extrapolation_bc = 0.0)

Calculate the combined bound-free cross section of an ion species using data from the Opacity
Project. The result includes contributions from all (provided) energy states, weighted by the
Boltzmann distribution, and includes the LTE correction for for stimulated emission.

# Arguments
- `λ_vals`: an iterable collection of wavelengths (in Å) to compute the cross sections at
- `T_vals`: an iterable collection of temperatures (in K) to compute the cross sections at
- `species_name`: name of the ion species

# Keyword Arguments
- `convert_to_cm2`: When True, returns the weighted cross section in units of cm². Otherwise, the
  results have units of megabarnes. Default is False.
- `partition_func`: Specifies the partition function for the current species. This should be a
  callable that accepts Temperature as an argument. When this is `nothing`, it falls back to the
  default partition functions.
- `extrapolation_bc`: Specifies handling of extrapolation during linear interpolation (this is
  passed to Interpolation.LinearInterpolation.

# Explanation
In more mathematical rigor, this function basically evaluates the following equation:

``\overline{\sigma_{{\rm bf},s}}(\lambda,T) = \sum_{{\bf n}\in{\rm states}} \frac{ g_{s,{\bf n}}}{Z_s(T)} e^{-E_{s,{\bf n}}/(k_B T)}\left(1 - e^{-h c/(\lambda k_B T)}\right) \sigma_{{\rm bf}, s,{\bf n}}(\lambda)``

where the ``s`` subscript corresponds to a given species. In this equation:
- ``g_{s,{\bf n}}`` is the statistical weight of state ``{\bf n}`` for species ``s``.
- ``Z_s`` is the partition function for species ``s``.
- ``E_{s,{\bf n}}`` is the excitation potential of state ``{\bf n}`` for species ``s``. In other
  words, its the energy difference between state ``{\bf n}`` and the ground state.
- ``\sigma_{\lambda,{\rm bf}, s,{\bf n}}`` is the bound free cross-section (a.k.a. photo-ionization
  cross section) of state ``{\bf n}`` for species ``s``.

Under the assumption of LTE, the linear aborption coefficient for bound-free absorption is simply:
``\alpha_{\lambda,{\rm bf},s} = n_{s,{\rm tot}} \overline{\sigma_{{\rm bf},s}}(\lambda,T)``.

# Notes
In the future, it needs to be possible to pass an argument that adjusts .
"""
function weighted_bf_cross_section_TOPBase(λ_vals, T_vals, species_name;
                                           elec_conf_table = nothing,
                                           cross_sec_file::Union{AbstractString,Nothing} = nothing,
                                           convert_to_cm2::Bool = false, partition_func = nothing,
                                           extrapolation_bc=0.0)

    # make it possible to pass in a dict mapping dict names to the various energy levels?
    # maybe also allow the user to specify an ionization energy level

    @assert ndims(λ_vals) == 1
    @assert ndims(T_vals) == 1

    # determine output units:
    units_factor = convert_to_cm2 ? 1e-18 : 1.0; # 1 megabarn = 1e-18 cm²

    # precompute Temperature-dependent constant
    inv_partition_func_val = if isnothing(partition_func)
        1.0./partition_funcs[Species(species_name)].(T_vals)
    else
        1.0./partition_func.(T_vals)
    end

    β_Ryd = RydbergH_eV./(kboltz_eV .* T_vals)

    # convert λ_vals to photon energies
    photon_energies = (hplanck_eV * c_cgs * 1.0e8) ./ λ_vals ./ RydbergH_eV

    # prepare the output array where results will be accumulated
    weighted_average = zeros(eltype(photon_energies), (length(photon_energies),length(T_vals)))

    # prepare the last few items before entering the loop
    state_prop_dict,itr = _get_tabulated_data(species_name)

    i = 0
    for (species, state_id, table_photon_energy_ryd, table_cross_section_MBarn) in itr
        i+=1
        if length(table_photon_energy_ryd) == 0
            continue
        end

        # construct an interpolator from the subtable
        func = _interpolator_from_subtable(table_photon_energy_ryd, table_cross_section_MBarn,
                                           species_name, extrapolation_bc)

        # retrieve the excitation energy (energy relative to ground state) and statistical weight
        state_prop = state_prop_dict[state_id]
        statistical_weight = state_prop.statistical_weight
        energy_ryd = abs(state_prop.excitation_potential_ryd)

        # now iterate over Temperatures
        for j in 1:length(T_vals)

            # now compute the weighting of the current state under LTE
            weight = inv_partition_func_val[j] * statistical_weight * exp(-energy_ryd * β_Ryd[j])

            current_cross_section = (func.(photon_energies) .*
                                     (1.0 .- exp.(-photon_energies.*β_Ryd[j]))
                                     .* units_factor)
            view(weighted_average, :, j) .+= (current_cross_section .* weight)
        end
    end

    weighted_average
end



"""
    absorption_coef_bf_TOPBase(λ, T, ndens_species, species_name; kwargs...)

Computes the bound-free linear absorption coefficient, α, (in units of cm⁻¹) using data from the 
Opacity Project (downloaded from TOPBase) and assuming LTE.

# Arguments
- `λ`: a scalar wavelength or an iterable collection of wavelengths to compute α at. This should 
  have units of ångstroms.
- `T`: a scalar temperature or an iterable collection of temperatures to compute α at. This should 
  have units of Kelvin.
- `ndens_species`: The number density of the species at each temperature (specified with units of
  cm⁻³). If `T` is a scalar, this should also be a scalar while if `T` is an iterable collection,
  this should have the same shape as `T`.
- `species_name`: name of the ion species

# Keyword Arguments
This function accepts all the same keyword arguments as `weighted_bf_cross_section_TOPBase`, except
for `convert_to_cm2`.
"""
function absorption_coef_bf_TOPBase(λ, T, ndens_species, species_name;
                                    kwargs...)

    for kwarg in keys(kwargs)
        if kwarg == :convert_to_cm2
            throw(ArgumentError(
                "absorption_coef_bf_TOPBase doesn't accept the convert_to_cm2 keyword argument"
            ))
        end
    end

    # the performance page of the Julia manual documents a design pattern with multiple dispatch
    # that could probably be used to acheive better performance
    @assert (ndims(λ) <= 1) & (ndims(T) <= 1)
    λ_vals = ndims(λ) == 0 ? typeof(λ)[λ] : λ

    @assert size(T) == size(ndens_species)
    (T_vals, ndens_species_vals) = if ndims(T) == 0
        (typeof(T)[T],                          # 1 element array holding T
         typeof(ndens_species)[ndens_species])  # 1 element array holding ndens_species
    else
        (T, ndens_species)
    end

    cross_sections = weighted_bf_cross_section_TOPBase(λ_vals, T_vals, species_name;
                                                       convert_to_cm2 = true, kwargs...)

    α = Array{eltype(cross_sections),2}(undef, size(cross_sections))
    for (i,n) in enumerate(ndens_species_vals)
        view(α,:,i) .= n .* view(cross_sections,:,i)
    end

    if ndims(λ) == 0 && ndims(T) == 0
        α[1,1]
    elseif ndims(λ) == 0
        α[1,:]
    elseif ndims(T) == 0
        α[:,1]
    else
        α
    end
end
