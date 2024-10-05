struct Wavelengths{F} <: AbstractArray{F,1}
    wl_ranges::Vector{StepRangeLen{F}} # in cm, not Å
    # these are for efficient-ish iteration, but ideally they should be eliminated
    all_wls::Vector{F}
    all_freqs::Vector{F}

    """
    TODO
    # Keyword Arguments
    - `air_wavelengths` (default: `false`): Whether or not the input wavelengths are air wavelengths to
        be converted to vacuum wavelengths by Korg.  The conversion will not be exact, so that the
        wavelength range can internally be represented by an evenly-spaced range.  If the approximation
        error is greater than `wavelength_conversion_warn_threshold`, an error will be thrown. (To do
        wavelength conversions yourself, see [`air_to_vacuum`](@ref) and [`vacuum_to_air`](@ref).)
    - `wavelength_conversion_warn_threshold` (default: 1e-4): see `air_wavelengths`. (In Å.)
    """
    function Wavelengths(wl_ranges::AbstractVector;
                         air_wavelengths=false, wavelength_conversion_warn_threshold=1e-4)
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

        new{eltype(all_wls)}(wl_ranges, all_wls, all_freqs)
    end
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
# handle integer args
function Wavelengths(λ_start::Integer, λ_stop::Integer, λ_step=0.01; kwargs...)
    Wavelengths(Float64(λ_start), Float64(λ_stop), λ_step; kwargs...)
end
function Wavelengths(λ_start::Integer, λ_stop, λ_step=0.01; kwargs...)
    Wavelengths(Float64(λ_start), λ_stop, λ_step; kwargs...)
end
function Wavelengths(λ_start, λ_stop::Integer, λ_step=0.01; kwargs...)
    Wavelengths(λ_start, Float64(λ_stop), λ_step; kwargs...)
end
# easy mode: pass in a single start and stop
function Wavelengths(λ_start, λ_stop, λ_step=0.01; kwargs...)
    Wavelengths([range(; start=λ_start, stop=λ_stop, step=λ_step)]; kwargs...)
end

# implement the AbstractArray interface
# https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array
Base.IndexStyle(::Type{<:Wavelengths}) = IndexLinear()
Base.size(wl::Wavelengths) = (length(wl.all_wls),) # implicitely defines Base.length
Base.getindex(wl::Wavelengths, i) = wl.all_wls[i]

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
TODO
"""
eachwindow(wls::Wavelengths) = ((first(r), last(r)) for r in wls.wl_ranges)

"""
TODO
TODO consider reversing?
"""
eachfreq(wls::Wavelengths) = wls.all_freqs

# index of the first element greater than or equal to λ
function Base.Sort.searchsortedfirst(wls::Wavelengths, λ)
    if λ >= 1 # convert Å to cm
        λ *= 1e-8
    end
    range_ind = searchsortedfirst(last.(wls.wl_ranges), λ)
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
    range_ind = searchsortedlast(first.(wls.wl_ranges), λ)
    if range_ind == 0
        return 0
    end
    ind = searchsortedlast(wls.wl_ranges[range_ind], λ)
    for ri in 1:range_ind-1
        ind += length(wls.wl_ranges[ri])
    end
    ind
end