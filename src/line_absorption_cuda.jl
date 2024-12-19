using SpecialFunctions: gamma
using ProgressMeter: @showprogress
using CUDA #TODO only import what you need

# TODO make it work with ForwardDiff.Dual
struct LineVals
    wl::Float64
    E_lower::Float64
    log_gf::Float64
    gamma_rad::Float64
    gamma_stark::Float64
    vdW::Tuple{Float64,Float64}
    ismolecular::Bool
    species_index::Int32

    function LineVals(line::Line, species_index::Int32)
        new(line.wl,
            line.E_lower,
            line.log_gf,
            line.gamma_rad,
            line.gamma_stark,
            if line.vdW isa Tuple
                line.vdW
            else
                (line.vdW, -1.0)
            end,
            ismolecule(line.species),
            species_index)
    end
end

"""
    line_absorption!(α, linelist, λs, temp, nₑ, n_densities, partition_fns, ξ
                   ; α_cntm=nothing, cutoff_threshold=1e-3, window_size=20.0*1e-8)

Calculate the opacity coefficient, α, in units of cm^-1 from all lines in `linelist`, at wavelengths
`λs` [cm^-1].

other arguments:

  - `temp` the temerature in K (as a vector, for multiple layers, if you like)
  - `n_densities`, a Dict mapping species to absolute number density in cm^-3 (as a vector, if temp is
    a vector).
  - `partition_fns`, a Dict containing the partition function of each species
  - `ξ` is the microturbulent velocity in cm/s (n.b. NOT km/s)
  - `α_cntm` is as a callable returning the continuum opacity as a function of wavelength. The window
    within which a line is calculated will extend to the wavelength at which the Lorentz wings or
    Doppler core of the line are at `cutoff_threshold * α_cntm[line.wl]`, whichever is greater.

# Keyword Arguments

  - `cuttoff_threshold` (default: 3e-4): see `α_cntm`
  - `verbose` (default: false): NOT ALLOWED TODO
"""
function line_absorption_cuda!(α, linelist, λs::Wavelengths, temps, nₑ, n_densities,
                               partition_fns, ξ, α_cntm; cutoff_threshold=3e-4, verbose=false)
    @assert verbose == false
    # allocate abs coef result array on the device
    α_d = CUDA.fill(0.0, size(α))
    line_absorption_cuda_helper!(α_d, linelist, λs, temps, nₑ, n_densities, partition_fns, ξ,
                                 α_cntm, cutoff_threshold)
    α .+= Array(α_d)
end

function line_absorption_cuda_helper!(α, linelist, λs::Wavelengths, temps, nₑ, n_densities,
                                      partition_fns, ξ, α_cntm, cutoff_threshold=3e-4,
                                      n_gpu_blocks=100)
    if length(linelist) == 0
        return zeros(length(λs))
    end

    each_species = unique((l.species for l in linelist))
    if species"H I" in each_species
        @error "Atomic hydrogen should not be in the linelist. Korg has built-in hydrogen lines."
    end
    species_index_dict = Dict(spec => Int32(i) for (i, spec) in enumerate(each_species))
    n_div_Z_cpu = zeros(eltype(α), (length(temps), length(each_species)))
    for (i, spec) in enumerate(each_species)
        n_div_Z_cpu[:, i] .= n_densities[spec] ./ partition_fns[spec].(log.(temps))
    end
    n_div_Z = CuArray(n_div_Z_cpu)
    mass_per_line_d = CuArray([get_mass(species) for species in each_species])

    # Each element of the arrays corresponds to an atmospheric layer, same at the "temps" array and 
    # the values in "number_densities"

    # these are line-independent
    temps_d = CuArray(temps)
    nₑ_d = CuArray(nₑ)
    n_H_I_d = CuArray(n_densities[species"H_I"])
    β = CuVector(@. 1 / (kboltz_eV * temps))
    λs_d = Wavelengths{eltype(λs.all_wls),
                       CuVector{eltype(λs.all_wls)},
                       eltype(λs.wl_ranges),
                       CuVector{eltype(λs.wl_ranges)}}(CuVector(λs.wl_ranges),
                                                       CuVector(λs.all_wls),
                                                       CuVector(λs.all_freqs))

    # these are line-dependent
    γ = CuArray{eltype(α)}(undef, n_gpu_blocks, length(temps))
    σ = CuArray{eltype(α)}(undef, n_gpu_blocks, length(temps))
    amplitude = CuArray{eltype(α)}(undef, n_gpu_blocks, length(temps))

    # convert the α_cntm interpolators to a matrix of coefficients on the coarse wavelength grid 
    # (i.e., don't interpolate)
    α_cntm_cpu = if eltype(α_cntm) <: Interpolations.AbstractInterpolation
        alphas = map(α_cntm) do itp
            itp.itp.coefs
        end
        vcat((a' for a in alphas)...)
    else
        # this is the 5000 Å only case
        reshape([a(0) for a in α_cntm], :, 1)
    end
    α_cntm_d = CuArray(α_cntm_cpu)
    coarse_λs_cpu = if eltype(α_cntm) <: Interpolations.AbstractInterpolation
        α_cntm[1].itp.knots[1].all_wls # TODO don't rely on the internals of Wavelengths?
    else
        [5000.0]
    end
    coarse_λs_d = CuArray(coarse_λs_cpu)
    all_line_vals_d = CuVector([LineVals(line, species_index_dict[line.species])
                                for (i, line) in enumerate(linelist)])
    warp_size = warpsize(device()) # need the device() call because this is called on CPU

    CUDA.@sync begin
        @cuda threads=warp_size blocks=n_gpu_blocks line_absorption_cuda_kernel!(α, all_line_vals_d,
                                                                                 σ,
                                                                                 λs_d, temps_d, β,
                                                                                 ξ, γ,
                                                                                 nₑ_d, n_H_I_d,
                                                                                 n_div_Z,
                                                                                 mass_per_line_d,
                                                                                 amplitude,
                                                                                 α_cntm_d,
                                                                                 coarse_λs_d,
                                                                                 cutoff_threshold)
    end
end

function line_absorption_cuda_kernel!(α, all_line_vals, σ, λs_d, temps_d, β, ξ, γ, nₑ_d, n_H_I_d,
                                      n_div_Z, mass_per_line_d, amplitude, α_cntm_d, coarse_λs_d,
                                      cutoff_threshold)
    for line_index in blockIdx().x:gridDim().x:length(all_line_vals)
        line_vals = all_line_vals[line_index]
        spec_index = line_vals.species_index
        process_line_kernel!(α, σ, λs_d, line_vals, temps_d, β, ξ, γ, nₑ_d, n_H_I_d, n_div_Z,
                             mass_per_line_d, amplitude, spec_index, α_cntm_d, coarse_λs_d,
                             cutoff_threshold)
    end
end

"""
This must be launched with threads equal to the warp size.
"""
function process_line_kernel!(α, σ, λs_d, line, temps_d, β, ξ, γ, nₑ_d, n_H_I_d,
                              n_div_Z, mass_per_line_d, amplitude, spec_index, α_cntm_d,
                              coarse_λs_d, cutoff_threshold)
    m = mass_per_line_d[spec_index]
    doppler_line_window = 0.0
    lorentz_line_window = 0.0
    blk_idx = blockIdx().x
    for idx in threadIdx().x:blockDim().x:size(σ, 2)
        σ[blk_idx, idx] = doppler_width(line.wl, temps_d[idx], m, ξ)

        # sum up the damping parameters.  These are FWHM (γ is usually the Lorentz HWHM) values in 
        # angular, not cyclical frequency (ω, not ν).
        Γ = line.gamma_rad
        if !line.ismolecular
            Γ += nₑ_d[idx] * scaled_stark(line.gamma_stark, temps_d[idx])
            Γ += n_H_I_d[idx] * scaled_vdW(line.vdW, m, temps_d[idx])
        end

        # calculate the lorentz broadening parameter in wavelength. Doing this involves an 
        # implicit aproximation that λ(ν) is linear over the line window.
        # the factor of λ²/c is |dλ/dν|, the factor of 1/2π is for angular vs cyclical freqency,
        # and the last factor of 1/2 is for FWHM vs HWHM
        γ[blk_idx, idx] = Γ * line.wl^2 / (c_cgs * 4π)

        E_upper = line.E_lower + c_cgs * hplanck_eV / line.wl
        levels_factor = exp(-β[idx] * line.E_lower) - exp(-β[idx] * E_upper)

        #total wl-integrated absorption coefficient
        amplitude[blk_idx, idx] = 10.0^line.log_gf * sigma_line(line.wl) *
                                  levels_factor * n_div_Z[idx, spec_index]

        #TODO don't do the searchsortedfirst at every layer?
        local_α_cntm = α_cntm_d[idx, searchsortedfirst(coarse_λs_d, line.wl)]
        ρ_crit = local_α_cntm * cutoff_threshold / amplitude[blk_idx, idx]

        inverse_gaussian_densities = inverse_gaussian_density(ρ_crit, σ[blk_idx, idx])
        inverse_lorentz_densities = inverse_lorentz_density(ρ_crit, γ[blk_idx, idx])

        doppler_line_window = max(doppler_line_window, inverse_gaussian_densities)
        lorentz_line_window = max(lorentz_line_window, inverse_lorentz_densities)
    end

    doppler_line_window = warp_reduce_max(doppler_line_window)
    lorentz_line_window = warp_reduce_max(lorentz_line_window)
    window_size = sqrt(lorentz_line_window^2 + doppler_line_window^2)

    lb = searchsortedfirst(λs_d, line.wl - window_size)
    ub = searchsortedlast(λs_d, line.wl + window_size)
    if lb > ub # TODO test performance
        return
    end

    for wl_index in lb:ub, thread_index in threadIdx().x:blockDim().x:size(σ, 2)
        # TODO DON'T BROADCAST
        Δα = line_profile.(line.wl, σ[blk_idx, thread_index], γ[blk_idx, thread_index],
                           amplitude[blk_idx, thread_index], λs_d[wl_index])
        ptr = pointer(α, thread_index + (wl_index - 1) * size(α, 1))
        CUDA.atomic_add!(ptr, Δα)
    end
    return
end

"""
This must be launched with threads equal to the warp size.
"""
function warp_reduce_max(value)
    offset = warpsize() ÷ 2 # it's not warpsize(device()) because this in part of a CUDA kernel
    while offset > 0
        value = max(value, CUDA.shfl_xor_sync(CUDA.FULL_MASK, value, offset))
        offset ÷= 2
    end
    return value
end