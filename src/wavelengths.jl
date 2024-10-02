struct Wavelengths{F} where R<:AbstractRange
    wl_ranges::Vector{R} # in cm, not Å
    # these are for efficient-ish iteration, but ideally they should be eliminated
    all_wls::Vector{F}
    all_freqs::Vector{F}

    """
    TODO
    """
    function Wavelengths(wl_ranges::AbstractVector{R}; air_wavelengths=false) where R
        # this could be more efficient
        all_λs = vcat(wl_ranges...)
        if !issorted(all_λs) #TODO test
            throw(ArgumentError("wl_ranges must be sorted and non-overlapping"))
        end

        #TODO handle air_wavelengths
        # Convert air to vacuum wavelenths if necessary.
        if air_wavelengths
            wl_ranges = map(wl_ranges) do wls
                λ_start, λ_stop, λ_step = first(wls), last(wls), step(wls)
                len = Int(round((λ_stop - λ_start) / λ_step)) + 1
                vac_start, vac_stop = air_to_vacuum.((λ_start, λ_stop))
                vac_step = (vac_stop - vac_start) / (len - 1)
                wls = StepRangeLen(vac_start, vac_step, len)
                max_diff = maximum(abs.(wls .- air_to_vacuum.(λ_start:λ_step:λ_stop)))
                if max_diff > wavelength_conversion_warn_threshold
                    throw(ArgumentError("A linear air wavelength range can't be approximated exactly with a"
                                        *
                                        "linear vacuum wavelength range. This solution differs by up to " *
                                        "$max_diff Å.  Adjust wavelength_conversion_warn_threshold if you " *
                                        "want to suppress this error."))
                end
                wls
            end
        end

        # TODO auto handle units?
        wl_ranges .*= 1e-8

        # precompute all wavelengths and frequencies
        all_wls = vcat(wl_ranges...)
        all_freqs = reverse(Korg.c_cgs ./ all_wls)

        new{R}(wl_ranges, all_wls, all_freqs)
    end
end
Wavelengths(wls::Wavelengths; kwargs...) = Wavelengths(wls.wl_ranges; kwargs...)
Wavelengths(wls::R; kwargs...) where R<:AbstractRange = Wavelengths([wls]; kwargs...)
function Wavelengths(wls::AbstractVector{<:Real}; tolerance=1e-6, kwargs...)
    min_step, max_step = extrema(diff(wls))
    if max_step - min_step > tolerance
        throw(ArgumentError("wavelengths are not linearly spaced to within $tolerance."))
    end
    Wavelengths([range(first(v), last(v); length=length(v))]; kwargs...)
end
function Wavelengths(λ_start, λ_stop, λ_step=0.01; kwargs...)
    Wavelengths([StepRangeLen(λ_start, λ_step, Int(round((λ_stop - λ_start) / λ_step)) + 1)];
                kwargs...)
end

"""
TODO
"""
eachrange(wls::Wavelengths) = wls.wl_range

"""
TODO
"""
eachwl(wl::Wavelengths) = wls.all_wls

"""
TODO
"""
eachfreq(wl::Wavelengths) = wls.all_freqs