"""
Functions for fitting to data.

!!! warning
    This submodule is in alpha. Do not use it for science.
"""
module Fit
using ..Korg, ProgressMeter, ForwardDiff, LineSearches, Optim
using SparseArrays: spzeros
using Statistics: mean
using Interpolations: LinearInterpolation
using InteractiveUtils: @code_warntype

# TODO specify versions in Project.toml
# TODO derivatives are zero at grid points (perturbing in one direction would be bad???)

"""
    compute_LSF_matrix(synth_wls, obs_wls, R; window_size=3)

Construct a sparse matrix, which when multiplied with a flux vector defined over wavelenths 
`synth_wls`, applies a gaussian line spead function (LSF) and resamples to the wavelenths `obswls`.
The LSF has a constant spectral resolution, ``R = \\lambda/\\Delta\\lambda``, where 
``\\Delta\\lambda`` is the LSF FWHM.  The `window_size` argument specifies how far out to extend
the convolution kernel in standard deviations.

For the best match to data, your wavelength range should extend a couple ``\\Delta\\lambda`` outside 
the region you are going to compare.

[`Korg.constant_R_LSF`](@ref) can apply an LSF to a single flux vector efficiently. This function is
relatively slow, but one the LSF matrix is constructed, convolving spectra to observational 
resolution via multiplication is fast.
"""
function compute_LSF_matrix(synth_wls, obs_wls, R; window_size=4)
    if !(first(synth_wls) < first(obs_wls) < last(obs_wls) < last(synth_wls))
        @warn raw"Synthesis wavelenths are not superset of observation wavelenths."
    end
    convM = spzeros((length(obs_wls), length(synth_wls)))
    lb, ub = 1,1 #initialize window bounds
    @showprogress  "Constructing LSF matrix" for i in 1:length(obs_wls)
        λ0 = obs_wls[i]
        σ = λ0 / R / (2sqrt(2log(2))) # convert Δλ = λ0/R (FWHM) to sigma
        lb, ub = Korg.move_bounds(synth_wls, lb, ub, λ0, window_size*σ)
        ϕ = Korg.normal_pdf.(synth_wls[lb:ub] .- λ0, σ) * step(synth_wls)
        @. convM[i, lb:ub] += ϕ
    end
    convM 
end

"""
Rescale each parameter so that it lives on (-∞, ∞) instead of the finite range of the atmosphere 
grid.
"""
function scale(p; lower = first.(Korg.get_atmosphere_archive()[1]),
                  upper = last.(Korg.get_atmosphere_archive()[1]))
    p = (p - lower) ./ (upper - lower)
    tan.(π .* (p.-0.5))
end

"""
Unscale each parameter so that it lives on the finite range of the atmosphere grid insteaf of 
(-∞, ∞).
"""
function unscale(p; lower = first.(Korg.get_atmosphere_archive()[1]),
                    upper = last.(Korg.get_atmosphere_archive()[1]))
    p = atan.(p)./π .+ 0.5
    p.*(upper - lower) .+ lower
end

"""
Synthesize a spectrum, returning the flux, with LSF applied, resampled, and rectified.  This is 
used by fitting routines. See [`Korg.synthesize`](@ref) to synthesize spectra as a Korg user.
"""
function synth(synthesis_wls, obs_wls, scaled_p::Vector{R}, linelist, LSF_matrix; line_buffer=10
              ) :: Vector{R} where R <: Real
    p = unscale(scaled_p)
    atm = Korg.interpolate_marcs(p...)
    alpha_els = ["Ne", "Mg", "Si", "S", "Ar", "Ca", "Ti"]
    A_X = Korg.format_A_X(p[3], Dict(["C"=>p[5]+p[3], (alpha_els .=> p[4]+p[3])...]))
    F::Vector{R} = try
        Korg.synthesize(atm, linelist, A_X, synthesis_wls; vmic=1.25, line_buffer=line_buffer).flux
    catch e
        ones(sum(length.(synthesis_wls)))
    end
    #rF = apply_rotation(F, wls, 10) #TODO
    dF = LSF_matrix * F
end

"""
TODO rotation, vmic, abundance

!!! warning
    This function is in alpha. Do not use it for science.
"""
function find_best_params_globally(obs_wls, obs_flux, obs_err, linelist, p0, synthesis_wls, LSF_matrix
                                  ; rectify=data_safe_rectify)
    obs_flux, obs_err, obs_wls, LSF_matrix = mask_out_nans(obs_flux, obs_err, obs_wls, LSF_matrix)
    rect_data = rectify(obs_flux, obs_err, obs_wls)
    rect_err = obs_err .* rect_data ./ obs_flux # divide out continuum
    chi2 = let obs_wls=obs_wls, data=rect_data, obs_err=obs_err, rect_err=rect_err, 
        synthesis_wls=synthesis_wls, LSF_matrix=LSF_matrix, linelist=linelist
        scaled_p -> begin
            flux = synth(synthesis_wls, obs_wls, scaled_p, linelist, LSF_matrix)
            flux = rectify(flux, obs_err, obs_wls)
            sum(((flux .- data)./rect_err).^2)
        end
    end 
    optimize(chi2, scale(p0), BFGS(linesearch=LineSearches.BackTracking()),  
            Optim.Options(x_tol=1e-5, time_limit=10_000, store_trace=true, 
            extended_trace=true), autodiff=:forward)
end

"""
TODO

using, e.g. iron lines
"""
function find_best_params_windowed(obs_wls, obs_flux, obs_err, windows, linelist, p0, synthesis_wls, 
                                   LSF_matrix; buffer=1.0)

    obs_flux, obs_err, obs_wls, LSF_matrix = mask_out_nans(obs_flux, obs_err, obs_wls, LSF_matrix)

    # calculate some synth ranges which span all windows, and the LSF submatrix that maps to them only
    windows = merge_bounds(windows, 2buffer)
    obs_wl_inds = map(windows) do (ll, ul)
        (findfirst(obs_wls .> ll), findfirst(obs_wls .> ul)-1)
    end
    obs_wl_mask, masked_obs_wls_range_inds, synth_wl_mask, multi_synth_wls = 
        calculate_multilocal_masks_and_ranges(obs_wl_inds, obs_wls, synthesis_wls, buffer)

    chi2 = let obs_wls=obs_wls[obs_wl_mask], data=obs_flux[obs_wl_mask], 
               obs_err=obs_err[obs_wl_mask], synthesis_wls=multi_synth_wls, 
               LSF_matrix=LSF_matrix[obs_wl_mask, synth_wl_mask], linelist=linelist
        scaled_p -> begin
            flux = synth(synthesis_wls, obs_wls, scaled_p, linelist, LSF_matrix)
            cntm = synth(synthesis_wls, obs_wls, scaled_p, [], LSF_matrix)
            flux ./= cntm
            sum(((flux .- data)./obs_err).^2)
        end
    end 
    
    res = optimize(chi2, scale(p0), BFGS(linesearch=LineSearches.BackTracking()),  
                   Optim.Options(x_tol=1e-5, time_limit=10_000, store_trace=true, 
                   extended_trace=true), autodiff=:forward)

    flux = synth(multi_synth_wls, obs_wls[obs_wl_mask], res.minimizer, linelist,
                 LSF_matrix[obs_wl_mask, synth_wl_mask])
    cntm = synth(multi_synth_wls, obs_wls[obs_wl_mask], res.minimizer, [],
                 LSF_matrix[obs_wl_mask, synth_wl_mask])
    best_fit_flux = flux ./ cntm

    res, obs_wl_mask, best_fit_flux
end

function find_best_fit_abundace(obs_wls, obs_flux, obs_err, linelist, windows, Z_to_fit, A_X_0, atm,
                                synthesis_wls, LSF_matrix; buffer=1.0)
    obs_flux, obs_err, obs_wls, LSF_matrix = mask_out_nans(obs_flux, obs_err, obs_wls, LSF_matrix)

    # calculate some synth ranges which span all windows, and the LSF submatrix that maps to them only
    windows = merge_bounds(windows, 2buffer)
    obs_wl_inds = map(windows) do (ll, ul)
        (findfirst(obs_wls .> ll), findfirst(obs_wls .> ul)-1)
    end
    obs_wl_mask, masked_obs_wls_range_inds, synth_wl_mask, multi_synth_wls = 
        calculate_multilocal_masks_and_ranges(obs_wl_inds, obs_wls, synthesis_wls, buffer)

    ll = A_X_0[Z_to_fit] - 4
    ul = A_X_0[Z_to_fit] + 4


    chi2 = let atm=atm, Z_to_fit=Z_to_fit, obs_wls=obs_wls[obs_wl_mask], 
               data=obs_flux[obs_wl_mask], obs_err=obs_err[obs_wl_mask], 
               synthesis_wls=multi_synth_wls, A_X_0=A_X_0,
               LSF_matrix=LSF_matrix[obs_wl_mask, synth_wl_mask], linelist=linelist
        scaled_abund -> begin
            abund = unscale(scaled_abund[1], lower=ll, upper=ul)
            A_X = [A_X_0[1:Z_to_fit-1] ; abund; A_X_0[Z_to_fit+1:length(A_X_0)]]
            flux = Korg.synthesize(atm, linelist, A_X, synthesis_wls).flux
            cntm = Korg.synthesize(atm, [], A_X, synthesis_wls).flux

            flux = LSF_matrix * (flux ./ cntm)
            sum(((flux .- data)./obs_err).^2)
        end
    end
    
    res = optimize(chi2, [scale(A_X_0[Z_to_fit], lower=ll, upper=ul)], autodiff=:forward,
                   Optim.Options(time_limit=10_000, store_trace=true, extended_trace=true))

    A_X = copy(A_X_0)
    A_X[Z_to_fit] = unscale(res.minimizer[1];  lower=ll, upper=ul)
    flux = Korg.synthesize(atm, linelist, A_X, multi_synth_wls).flux
    cntm = Korg.synthesize(atm, [], A_X, multi_synth_wls).flux
    best_fit_flux = LSF_matrix[obs_wl_mask, synth_wl_mask] * (flux ./ cntm)

    res, A_X[Z_to_fit], obs_wl_mask, best_fit_flux
end

"""
TODO

synthesize spectrum in windows which are ± wl_buffer larger than the regions used for 
comparison.  This is important because the spectra get convolved with the LSF. 

TODO local_rectify signature: error or not

TODO synthesis wls must be a range. This isn't a bug problem because the full-spectrum synthesis is
only done once.

!!! warning
    This function is in alpha. Do not use it for science.
"""
function find_best_params_multilocally(obs_wls, obs_flux, obs_err, linelist, p0, synthesis_wls, 
                                       LSF_matrix; verbose=true, global_rectify=data_safe_rectify, 
                                       local_rectify=simple_rectify, wl_chunk_size=10, wl_buffer=3,
                                       line_buffer=3, prerectified=false)
    obs_flux, obs_err, obs_wls, LSF_matrix = mask_out_nans(obs_flux, obs_err, obs_wls, LSF_matrix)
    J0 = global_jacobian(obs_wls, obs_err, linelist, p0, synthesis_wls, LSF_matrix)

    # contains *indices* into obs_wls for each wl range
    obs_bounds_inds = find_sensitive_ranges(J0, obs_err, obs_wls, wl_chunk_size)
    obs_bounds_inds = merge_bounds(obs_bounds_inds, 2wl_buffer / step(synthesis_wls))
    # The various masks and ranges we need to work with the subspectra. See docstring for details.
    obs_wl_mask, masked_obs_wls_range_inds, synth_wl_mask, multi_synth_wls = 
        calculate_multilocal_masks_and_ranges(obs_bounds_inds, obs_wls, synthesis_wls, wl_buffer)
    if verbose
        println("Using subspectra: ", ("$(round(obs_wls[lb], digits=1))–$(round(obs_wls[ub], digits=1)) "
                                       for (lb, ub) in obs_bounds_inds)...)
    end

    # rectify subspectra of input data
    rectified_local_data, rectified_local_err = if prerectified
        obs_flux[obs_wl_mask], obs_err[obs_wl_mask]
    else
        copy(obs_flux[obs_wl_mask]), copy(obs_err[obs_wl_mask])
        for r in masked_obs_wls_range_inds
            rF = local_rectify(rectified_local_data[r], obs_wls[obs_wl_mask][r], obs_err[obs_wl_mask][r])
            rectified_local_err[r] ./= rectified_local_data[r] ./ rF # divide cntm out of err
            rectified_local_data[r] .= rF
        end
    end

    function chi2(scaled_p) 
        flux = synth(multi_synth_wls, obs_wls[obs_wl_mask], scaled_p, linelist, 
                     LSF_matrix[obs_wl_mask, synth_wl_mask]; line_buffer=line_buffer)
        if prerectified
            cntm = synth(multi_synth_wls, obs_wls[obs_wl_mask], scaled_p, [], 
                         LSF_matrix[obs_wl_mask, synth_wl_mask]; line_buffer=line_buffer)
            flux ./= cntm
        else
            for r in masked_obs_wls_range_inds
                flux[r] .= local_rectify(flux[r], obs_wls[obs_wl_mask][r], obs_err[obs_wl_mask][r])
            end
        end
        sum(((flux .- rectified_local_data)./rectified_local_err).^2)
    end 
    result = optimize(chi2, scale(p0), BFGS(linesearch=LineSearches.BackTracking()),  
                            Optim.Options(time_limit=3600, x_tol=1e-5, store_trace=true, 
                            extended_trace=true), autodiff=:forward)

    (result, chi2, multi_synth_wls, LSF_matrix[obs_wl_mask, synth_wl_mask], obs_wls[obs_wl_mask], 
    rectified_local_data, rectified_local_err, masked_obs_wls_range_inds)
end

"""
Sort a vector of lower-bound, uppoer-bound pairs and merge overlapping ranges.  Used by 
find_best_params_multilocally.
"""
function merge_bounds(bounds, merge_distance)
    bounds = sort(bounds, by=first)
    new_bounds = [bounds[1]]
    for i in 2:length(bounds)
        # if these bounds are within merge_distance of the previous, extend the previous, 
        # otherwise add them to the list
        if bounds[i][1] < new_bounds[end][2] + merge_distance 
            new_bounds[end] = (new_bounds[end][1], max(bounds[i][2], new_bounds[end][2]))
        else
            push!(new_bounds, bounds[i])
        end
    end
    new_bounds
end

"""
Divide `flux` by its mean.  This is a good way to "rectify" small wavelength regions.
"""
function simple_rectify(flux, wls, err)
    iv = err.^-2
    meanF = sum(flux .* iv) / sum(iv)
    flux / meanF
end

"""
    data_safe_rectify(flux::Vector{R}, err, wls; bandwidth=50, wl_step=0)

Rectify the spectrum with flux vector `flux` and wavelengths `wls` by dividing out a moving
`q`-quantile with window size `bandwidth`.  `wl_step` controls the size of the grid that the moving 
quantile is calculated on and interpolated from. By default, this calculation is exact, but Setting
`wl_step` to a non-zero value results in calculating the moving mean every `wl_step` Å and 
interpolating to get the "continuum".

This method is not a true rectification in the sense that it does not set the continuum to 1, but
it removes the continuum from a spectrum without introducing any bias. It does not work well with 
linelists lacking many lines.
"""
function data_safe_rectify(flux::Vector{R}, err, wls; bandwidth=50, wl_step=0
                          ) :: Vector{R} where R <: Real
    inv_var = err.^-2
    #construct a range of integer indices into wls corresponding to roughly wl_step-sized steps
    if wl_step == 0
        inds = eachindex(wls)
    else
        inds = 1 : max(1, Int(floor(wl_step/step(wls)))) : length(wls)
    end
    lb = 1
    ub = 1
    moving_mean = map(wls[inds]) do λ
        lb, ub = Korg.move_bounds(wls, lb, ub, λ, bandwidth)
        sum(flux[lb:ub].*inv_var[lb:ub]) / sum(inv_var[lb:ub])
    end
    if wl_step == 0
        flux ./ moving_mean
    else
        itp = LinearInterpolation(wls[inds], moving_quantile, extrapolation_bc=Flat())
        flux ./ itp.(wls)
    end
end


"""
Given a Jacobian, `J`, whose data index is over wavelengths, `obs_wls`, find the 
`chunks_per_param` wavelength chunks of length `chunk_size` which are most sensitive to each 
parameter.
"""
function find_sensitive_ranges(J, err, obs_wls, chunk_size; chunks_per_param=3)
    # chucks overlap by 50% with their neighbors
    wl_half_chunk_edges = first(obs_wls) : chunk_size/2 : last(obs_wls)

    bounds = Vector{Tuple{Int, Int}}()
    lb, ub = 1, 1
    for ∂F_∂x in eachcol(J)
        local_grad2 = []
        for wl_ind in eachindex(wl_half_chunk_edges[1:end-1])
            wl_lb = wl_half_chunk_edges[wl_ind]
            wl_ub = wl_half_chunk_edges[wl_ind+1]
            lb, ub = Korg.move_bounds(obs_wls, lb, ub, (wl_lb+wl_ub)/2, chunk_size)

            # if there are gaps in the observed spectra, there may be no data in the wl range
            if lb < ub 
                push!(local_grad2, (lb, ub, sum(@. (∂F_∂x[lb:ub]/err[lb:ub])^2)))
            end
        end
        for i in partialsortperm(last.(local_grad2), 1:chunks_per_param, rev=true)
            lb, ub, m = local_grad2[i]
            #println(obs_wls[lb], "-", obs_wls[ub], " ", lb, " ", ub, " ", m)
            push!(bounds, (lb, ub))
        end
    end
    bounds
end

"""
    calculate_multilocal_masks_and_ranges(obs_bounds_inds, obs_wls, synthesis_wls)

Given a vector of target synthesis ranges in the observbed spectrum, return the masks, etc required.

Arguments:
    - `obs_bounds_inds`: a vector of tuples of the form `(lower_bound_index, upper_bound_index)`
    - `obs_wls`: the wavelengths of the observed spectrum
    - `synthesis_wls`: the wavelengths of the synthesis spectrum
    - `wl_buffer`: the number of Å to add to each side of the synthesis range

Returns:
    - `obs_wl_mask`: a bitmask for `obs_wls` which selects the observed wavelengths
    - `masked_obs_wls_range_inds`: a vector of ranges which select each subspectrum in the masked 
       observed spectrum. 
    - `synthesis_wl_mask`: a bitmask for `synthesis_wls` which selects the synthesis wavelengths 
       needed to generated the masked observed spectrum.
    - `multi_synth_wls`: The vector of ranges to pass to `Korg.synthesize`.
"""
function calculate_multilocal_masks_and_ranges(obs_bounds_inds, obs_wls, synthesis_wls, wl_buffer)
    obs_wl_mask = zeros(Bool, length(obs_wls)) 
    for (lb, ub) in obs_bounds_inds
        obs_wl_mask[lb:ub] .= true
    end
    # these ranges each select a subsprectrum from obs_wls[obs_wl_mask]
    masked_obs_wls_range_inds = [(1 : obs_bounds_inds[1][2] - obs_bounds_inds[1][1] + 1)]
    for (lb, ub) in obs_bounds_inds[2:end]
        prev_ub = masked_obs_wls_range_inds[end][end]
        push!(masked_obs_wls_range_inds, prev_ub+1 : prev_ub+ub-lb+1)
    end
    # bitmask for synthesis_wls to isolate the subspectra
    synth_wl_mask = zeros(Bool, length(synthesis_wls)) 
    # multi_synth_wls is the vector of wavelength ranges that gets passed to synthesize
    multi_synth_wls = map(obs_bounds_inds) do (lb, ub)
        synth_wl_lb = findfirst(synthesis_wls .> obs_wls[lb] - wl_buffer)
        synth_wl_ub = findfirst(synthesis_wls .> obs_wls[ub] + wl_buffer) - 1
        synth_wl_mask[synth_wl_lb:synth_wl_ub] .= true
        synthesis_wls[synth_wl_lb:synth_wl_ub]
    end
    obs_wl_mask, masked_obs_wls_range_inds, synth_wl_mask, multi_synth_wls
end

"""
TODO
"""
function mask_out_nans(obs_flux, obs_err, obs_wls, LSF_matrix)
    mask = @. !isnan(obs_flux) && !isnan(obs_err)
    obs_flux[mask], obs_err[mask], obs_wls[mask], LSF_matrix[mask, :]
end

"""
TODO
"""
function global_jacobian(obs_wls, obs_err, linelist, p0, synthesis_wls, LSF_matrix;
                         verbose=true, global_rectify=data_safe_rectify)
    J0 = let synthesis_wls=synthesis_wls, obs_wls=obs_wls, linelist=linelist, LSF_matrix=LSF_matrix, 
             obs_err=obs_err
        ForwardDiff.jacobian(scale(p0)) do p
            flux = synth(synthesis_wls, obs_wls, p, linelist, LSF_matrix)
            global_rectify(flux, obs_err, obs_wls)
        end
    end
    @assert !any(isnan.(J0))
    J0
end

end # module