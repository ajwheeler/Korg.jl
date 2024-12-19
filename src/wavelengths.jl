import Adapt

struct Wavelengths{F,V<:AbstractVector{F},R<:AbstractRange{F},VR<:AbstractVector{R}} <:
       AbstractArray{F,1}
    wl_ranges::VR # in cm, not Å
    # precomputed arrays.  This is faster for indexing into Wavelengths, especially with a UnitRange
    all_wls::V
    all_freqs::V # this one should probably be dropped.
end
"""
Korg.Wavelengths(wl_params...; air_wavelengths=false, wavelength_conversion_warn_threshold=1e-4)

Construct a `Wavelengths` object which represents the (possibly non-contiguous) wavelengths for
which to compute a spectrum.  The wavelengths can be specified with an upper and lower bound, or
a vector of upper and lower bounds. For example,

Korg.Wavelengths(5000, 5500)
Korg.Wavelengths([(5000, 5500), (6000, 6500)])

# Keyword Arguments

  - `air_wavelengths` (default: `false`): Whether or not the input wavelengths are air wavelengths to
    be converted to vacuum wavelengths by Korg.  The conversion will not be exact, so that the
    wavelength range can internally be represented by an evenly-spaced range.  If the approximation
    error is greater than `wavelength_conversion_warn_threshold`, an error will be thrown. (To do
    wavelength conversions yourself, see [`air_to_vacuum`](@ref) and [`vacuum_to_air`](@ref).)
  - `wavelength_conversion_warn_threshold` (default: 1e-4): see `air_wavelengths`. (In Å.)
"""
function Wavelengths(wl_ranges::AbstractVector{R};
                     air_wavelengths=false,
                     wavelength_conversion_warn_threshold=1e-4) where R<:AbstractRange
    # if the first wavelength is > 1, assume it's in Å and convert to cm
    if first(first(wl_ranges)) >= 1
        wl_ranges = wl_ranges .* 1e-8
    end

    # this could be more efficient
    all_wls = vcat(wl_ranges...)
    if !issorted(all_wls)
        throw(ArgumentError("wl_ranges must be sorted and non-overlapping"))
    end

    # Convert air to vacuum wavelenths if necessary.
    if air_wavelengths
        wl_ranges = map(wl_ranges) do wls
            vac_start, vac_stop = air_to_vacuum.((first(wls), last(wls)))
            vac_wls = range(; start=vac_start, stop=vac_stop, length=length(wls))
            max_diff = maximum(abs.(vac_wls .- air_to_vacuum.(wls)))
            if max_diff > wavelength_conversion_warn_threshold
                throw(ArgumentError("A linear air wavelength range can't be approximated exactly with a"
                                    *
                                    "linear vacuum wavelength range. This solution differs by up to " *
                                    "$max_diff Å.  Adjust wavelength_conversion_warn_threshold if you " *
                                    "want to suppress this error."))
            end
            vac_wls
        end
    end

    # precompute all wavelengths and frequencies
    all_freqs = reverse(Korg.c_cgs ./ all_wls)

    Wavelengths{eltype(all_wls),
                typeof(all_wls),
                eltype(wl_ranges),
                typeof(wl_ranges)}(wl_ranges, all_wls, all_freqs)
end
function Wavelengths(wls::Wavelengths; air_wavelengths=false,)
    if air_wavelengths
        Wavelengths(wls.wl_ranges; air_wavelengths=true, kwargs...)
    else
        wls
    end
end
Wavelengths(wls::R; kwargs...) where R<:AbstractRange = Wavelengths([wls]; kwargs...)
function Wavelengths(wls::AbstractVector{<:Real}; tolerance=1e-6, kwargs...)
    if length(wls) == 0
        throw(ArgumentError("wavelengths must be non-empty"))
    elseif length(wls) == 1
        Wavelengths(wls[1]:1.0:wls[1]; kwargs...)
    else
        min_step, max_step = extrema(diff(wls))
        if max_step - min_step > tolerance
            throw(ArgumentError("wavelengths are not linearly spaced to within $tolerance."))
        end
        Wavelengths([range(first(wls), last(wls); length=length(wls))]; kwargs...)
    end
end
function Wavelengths(wls::AbstractVector, λ_step=(wls[1][1] < 1 ? 0.01e-8 : 0.01); kwargs...)
    if λ_step isa Integer
        λ_step = convert(Float64, λ_step)
    end
    Wavelengths([range(; start=λ_start, stop=λ_stop, step=λ_step) for (λ_start, λ_stop) in wls];
                kwargs...)
end
function Wavelengths(λ_start, λ_stop, args...; kwargs...)
    Wavelengths([(λ_start, λ_stop)], args...; kwargs...)
end
Wavelengths(wls::Tuple{<:Real,<:Real}; kwargs...) = Wavelengths([wls]; kwargs...)

# implement the AbstractArray interface
# https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array
Base.IndexStyle(::Type{<:Wavelengths}) = IndexLinear() #1-dimensional indexing is "native" 
Base.size(wl::Wavelengths) = (length(wl.all_wls),) # implicitely defines Base.length
Base.getindex(wl::Wavelengths, i) = wl.all_wls[i]

function Adapt.adapt_structure(to, wls::Wavelengths)
    wl_ranges = Adapt.adapt_structure(to, wls.wl_ranges)
    all_wls = Adapt.adapt_structure(to, wls.all_wls)
    all_freqs = Adapt.adapt_structure(to, wls.all_freqs)
    Wavelengths{eltype(all_wls),
                typeof(all_wls),
                eltype(wl_ranges),
                typeof(wl_ranges)}(wl_ranges, all_wls, all_freqs)
end

function Base.show(io::IO, wl::Wavelengths)
    print(io, "Korg.Wavelengths(")
    for r in wl.wl_ranges[1:end-1]
        print(io, Int(round(first(r * 1e8))), "—", Int(round(last(r * 1e8))), ", ")
    end
    println(io, Int(round(first(wl.wl_ranges[end] * 1e8))), "—",
            Int(round(last(wl.wl_ranges[end]) * 1e8)), ")")
end
Base.show(io::IO, ::MIME"text/plain", wl::Wavelengths) = show(io, wl) # REPL/notebook

Base.:(==)(wl1::Wavelengths, wl2::Wavelengths) = wl1.wl_ranges == wl2.wl_ranges
function Base.isapprox(wl1::Wavelengths, wl2::Wavelengths; kwargs...)
    isapprox(wl1.all_wls, wl2.all_wls; kwargs...)
end

"""
    eachwindows(wls::Wavelengths)

Returns an iterator over the wavelength ranges `(λ_low, λ_hi)` in `wls` (in cm, not Å).
"""
eachwindow(wls::Wavelengths) = ((first(r), last(r)) for r in wls.wl_ranges)

"""
    eachfreq(wls::Wavelengths)

Returns an array of the frequencies corresponding to the wavelengths in `wls` (in Hz). They are
sorted, i.e. in reverse order of the wavelengths.
"""
eachfreq(wls::Wavelengths) = wls.all_freqs

"""
    subspectrum_indices(wls::Wavelengths)

Returns a vector of Julia ranges, which can be used to index into the full spectrum to get the
sub-spectrum corresponding to each wavelength range in `wls`.
"""
function subspectrum_indices(wls::Wavelengths)
    # collect the indices corresponding to each wavelength range
    wl_lb_ind = 1 # the index into α of the lowest λ in the current wavelength range
    indices = []
    for λs in wls.wl_ranges
        wl_inds = wl_lb_ind:wl_lb_ind+length(λs)-1
        push!(indices, wl_inds)
        wl_lb_ind += length(λs)
    end
    indices
end

# index of the first element greater than or equal to λ
function Base.Sort.searchsortedfirst(wls::Wavelengths, λ)
    if λ >= 1 # convert Å to cm
        λ *= 1e-8
    end
    range_ind = searchsortedfirst(wls.wl_ranges, λ; by=last)
    if range_ind == length(wls.wl_ranges) + 1
        return length(wls) + 1
    end
    ind = searchsortedfirst(wls.wl_ranges[range_ind], λ)
    for ri in 1:range_ind-1
        ind += length(wls.wl_ranges[ri])
    end
    ind
end
# index of the last element less than or equal to λ
function Base.Sort.searchsortedlast(wls::Wavelengths, λ)
    if λ >= 1 # convert Å to cm
        λ *= 1e-8
    end
    range_ind = searchsortedlast(wls.wl_ranges, λ; by=first)
    if range_ind == 0
        return 0
    end
    ind = searchsortedlast(wls.wl_ranges[range_ind], λ)
    for ri in 1:range_ind-1
        ind += length(wls.wl_ranges[ri])
    end
    ind
end