"""
Functions for fitting to data.

!!! warning

    This submodule is in beta. It's API may change.
"""
module Fit
using Compat: @compat
@compat public fit_spectrum# , ews_to_abundances, ews_to_stellar_parameters

using ..Korg, LineSearches, Optim
using Interpolations: linear_interpolation, Line
using ForwardDiff, DiffResults
using Trapz
using Statistics: mean, std
using ProgressMeter
using NonlinearSolve
using SciMLBase: successful_retcode
using FiniteDiff

# used by scale and unscale for some parameters
function tan_scale(p, lower, upper)
    if !(lower <= p <= upper)
        throw(ArgumentError("p=$p is not in the range $lower to $upper"))
    end
    tan(π * (((p - lower) / (upper - lower)) - 0.5))
end
tan_unscale(p, lower, upper) = (atan(p) / π + 0.5) * (upper - lower) + lower

# these are the parameters which are scaled by tan_scale
const tan_scale_params = Dict("epsilon" => (0, 1),
                              "cntm_offset" => (-0.5, 0.5),
                              "cntm_slope" => (-0.1, 0.1),
                              # we can't get these directly from Korg.get_atmosphere_archive() because it will fail in the
                              # test environment, but they are simply the boundaries of the SDSS marcs grid used by
                              # Korg.interpolate_marcs.
                              "Teff" => (2800, 8000),
                              "logg" => (-0.5, 5.5),
                              "m_H" => (-5, 1),
                              # this allows all the atmospheres supported by the grid, but also many that are not.
                              # alpha will be clamped to the nearest supported value.
                              "alpha_H" => (-3.5, 2),
                              map(Korg.atomic_symbols) do el
                                  el => (-10, +4)
                              end...)

"""
Rescale each parameter so that it lives on (-∞, ∞).
"""
scale(params::Dict) = map(collect(params)) do (name, p)
    name => if name in keys(tan_scale_params)
        tan_scale(p, tan_scale_params[name]...)
    elseif name in ["vmic", "vsini"]
        tan_scale(sqrt(p), 0, sqrt(250))
    else
        @error "$name is not a parameter I know how to scale."
    end
end |> Dict

"""
Unscale each parameter so that it lives on the appropriate range instead of (-∞, ∞).
"""
function unscale(params::Dict)
    map(collect(params)) do (name, p)
        name => if name in keys(tan_scale_params)
            tan_unscale(p, tan_scale_params[name]...)
        elseif name in ["vmic", "vsini"]
            tan_unscale(p, 0, sqrt(250))^2
        else
            @error "$name is not a parameter I know how to unscale."
        end
    end |> Dict
end

"""
Synthesize a spectrum, returning the flux, with LSF applied, resampled, and rectified.  This is
an internal function oned by fitting routines.
See [`Korg.synthesize`](@ref) to synthesize spectra as a Korg user.
"""
function synthetic_spectrum(synthesis_wls, linelist, LSF_matrix, params, synthesis_kwargs)
    specified_abundances = Dict([p for p in pairs(params) if p.first in Korg.atomic_symbols])
    alpha_H = "alpha_H" in keys(params) ? params["alpha_H"] : params["m_H"]
    A_X::Vector{valtype(params)} = Korg.format_A_X(params["m_H"], alpha_H, specified_abundances;
                                                   solar_relative=true)

    # clamp_abundances clamps M_H, alpha_M, and C_M to be within the atm grid
    atm = Korg.interpolate_marcs(params["Teff"], params["logg"], A_X; clamp_abundances=true,
                                 perturb_at_grid_values=true)

    sol = Korg.synthesize(atm, linelist, A_X, synthesis_wls; vmic=params["vmic"], line_buffer=0,
                          electron_number_density_warn_threshold=Inf, synthesis_kwargs...)

    # apply cntm adjustments
    central_wavelength = (sol.wavelengths[begin] + sol.wavelengths[end]) / 2
    cntm_adjustment = 1 .- params["cntm_offset"] .-
                      params["cntm_slope"] * (sol.wavelengths .- central_wavelength)
    F = sol.flux ./ (sol.cntm .* cntm_adjustment)

    F = Korg.apply_rotation(F, synthesis_wls, params["vsini"], params["epsilon"])
    LSF_matrix * F
end

"""
Synthesize a spectrum, apply the LSF, and postprocess it, catching and potentially rethrowing
errors. This is used by [`fit_spectrum`](@ref).
"""
function postprocessed_synthetic_spectrum(synth_wls, linelist, LSF_matrix, params,
                                          synthesis_kwargs, obs_wls, windows, obs_flux, obs_err,
                                          postprocess, adjust_continuum)
    flux = try
        synthetic_spectrum(synth_wls, linelist, LSF_matrix, params, synthesis_kwargs)
    catch e
        if (e isa Korg.ChemicalEquilibriumError) || (e isa Korg.LazyMultilinearInterpError)
            # chemical equilibrium errors happen for a few unphysical model atmospheres

            # LazyMultilinearInterpError happens when logg is oob for the low-Z atmosphere grid

            # This is a nice huge chi2 value, but not too big.  It's what you get if
            # difference at each pixel in the (rectified) spectra is 1, which is
            # more-or-less an upper bound.
            return sum(1 ./ obs_err .^ 2)
        else
            rethrow(e)
        end
    end

    try
        postprocess(flux, obs_flux, obs_err)
    catch e
        println(stderr, "Error while calling postprocess")
        rethrow(e)
    end

    if adjust_continuum
        try
            linear_continuum_adjustment!(obs_wls, windows, flux, obs_flux, obs_err)
        catch e
            println(stderr, "Error while calling adjust_continuum")
            rethrow(e)
        end
    end

    flux
end

"""
Validate fitting parameters, and insert default values when needed. Used by [`fit_spectrum`](@ref).

these can be specified in either initial_guesses or fixed_params, but if they are not, these values
are inserted into fixed_params
"""
function validate_params(initial_guesses::AbstractDict, fixed_params::AbstractDict;
                         required_params=["Teff", "logg"],
                         default_params=Dict("m_H" => 0.0, "vsini" => 0.0, "vmic" => 1.0,
                                             "epsilon" => 0.6, "cntm_offset" => 0.0,
                                             "cntm_slope" => 0.0),
                         allowed_params=Set(["alpha_H"; required_params; keys(default_params)...;
                                             Korg.atomic_symbols]))
    # convert all parameter values to Float64
    initial_guesses = Dict(string(p[1]) => Float64(p[2]) for p in pairs(initial_guesses))
    fixed_params = Dict(string(p[1]) => Float64(p[2]) for p in pairs(fixed_params))

    # check that all required params are specified
    all_params = keys(initial_guesses) ∪ keys(fixed_params)
    for param in required_params
        if !(param in all_params)
            throw(ArgumentError("Must specify $param in either initial_guesses or fixed_params. (Did you get the capitalization right?)"))
        end
    end

    # check that all the params are recognized
    unknown_params = filter!(all_params) do param
        param ∉ allowed_params
    end
    if length(unknown_params) > 0
        throw(ArgumentError("These parameters are not recognized: $(unknown_params)"))
    end

    # filter out those that the user specified, and fill in the rest
    default_params = filter(collect(pairs(default_params))) do (k, v)
        !(k in keys(initial_guesses)) && !(k in keys(fixed_params))
    end |> Dict
    fixed_params = merge(default_params, fixed_params)

    # check that no params are both fixed and initial guesses
    let keys_in_both = collect(keys(initial_guesses) ∩ keys(fixed_params))
        if length(keys_in_both) > 0
            throw(ArgumentError("These parameters: $(keys_in_both) are specified as both initial guesses and fixed params."))
        end
    end

    if "cntm_offset" in keys(initial_guesses) || "cntm_slope" in keys(initial_guesses)
        @warn "Instead of using the `\"cntm_offset\"`` and `\"cntm_slope\"` parameters, it's now" *
              " recommended to pass `adjust_continuum=true` to Korg.Fit.fit_spectrum. These parameters " *
              "may be removed in a future version of Korg."
    end

    initial_guesses, fixed_params
end

# make it possible to use dicts instead of NamedTuples for the python people
function validate_params(initial_guesses::AbstractDict, fixed_params::NamedTuple; kwargs...)
    validate_params(initial_guesses, _namedtuple_to_dict(fixed_params); kwargs...)
end
function validate_params(initial_guesses::NamedTuple, fixed_params=AbstractDict{String,Float64}();
                         kwargs...)
    validate_params(_namedtuple_to_dict(initial_guesses), fixed_params; kwargs...)
end

function _namedtuple_to_dict(nt::NamedTuple)
    Dict{String,Float64}([string(p[1]) => Float64(p[2]) for p in pairs(nt)])
end

"""
    fit_spectrum(obs_wls, obs_flux, obs_err, linelist, initial_guesses, fixed_params; kwargs...)

Find the parameters and abundances that best match a rectified observed spectrum.

# Arguments:

  - `obs_wls`: the wavelengths of the observed spectrum in Å.  These must be vacuum wavelengths.
  - `obs_flux`: the rectified flux of the observed spectrum
  - `obs_err`: uncertainty in `flux`
  - `linelist`: a linelist to use for the synthesis
  - `initial_guesses`: a NamedTuple specifying initial guesses for the parameters to be fit.  See
    "Specifying parameters" below.
  - `fixed_params`: a NamedTuple specifying parameters to be held fixed. See "Specifying parameters"
    below.

`initial_guesses` and `fixed_params` can also be specified as Dicts instead of NamedTuples, which is
more convenient when calling Korg from python.

# Specifying parameters

Parameters are specified as named tuples or dictionaries. Named tuples look like this:
`(Teff=5000, logg=4.5, m_H=0.0)`.  Single-element named tuples require a semicolon: `(; Teff=5000)`.

### Required parameters

`Teff` and `logg` *must* be specified in either `initial_guesses` or `fixed_params`.

### Optional Parameters

These can be specified in either `initial_guesses` or `fixed_params`, but if they are not default
values are used.

  - `m_H`: the metallicity of the star, in dex. Default: `0.0`
  - `alpha_H`: the alpha enhancement of the star, in dex. Default: `m_H`.  Note that, because of the
    parameter range supported by [`Korg.interpolate_marcs`](@ref), only values within ±1 of `m_H`
    are supported.
  - `vmic`: the microturbulence velocity, in km/s. Default: `1.0`
  - `vsini`: the projected rotational velocity of the star, in km/s. Default: `0.0`.
    See [`Korg.apply_rotation`](@ref) for details.
  - `epsilon`: the linear limb-darkening coefficient. Default: `0.6`. Used for applying rotational
    broadening only.  See [`Korg.apply_rotation`](@ref) for details.
  - Individual elements, e.g. `Na`, specify the solar-relative ([X/H]) abundance of that element.

# Keyword arguments

  - `R`, the resolution of the observed spectrum. This is required.  It can be specified as a
    function of wavelength, in which case it will be evaluated at the observed wavelengths.
  - `windows` is a vector of wavelength pairs, each of which specifies a wavelength
    "window" to synthesize and contribute to the total χ². If not specified, the entire spectrum is
    used. Overlapping windows are automatically merged.
  - `adjust_continuum` (default: `false`) if true, adjust the continuum with the best-fit linear
    correction within each window, minimizing the chi-squared between data and model at every step
    of the optimization.
  - `wl_buffer` is the number of Å to add to each side of the synthesis range for each window.
  - `time_limit` is the maximum number of seconds to spend in the optimizer. (default: `10_000`).
    The optimizer will only checks against the time limit after each step, so the actual wall time
    may exceed this limit.
  - `precision` specifies the tolerance for the solver to accept a solution. The solver operates on
    transformed parameters, so `precision` doesn't translate straightforwardly to Teff, logg, etc, but
    the default value, `1e-4`, provides a theoretical worst-case tolerance of about 0.15 K in `Teff`,
    0.0002 in `logg`, 0.0001 in `m_H`, and 0.0004 in detailed abundances. In practice the precision
    achieved by the optimizer is about 10x bigger than this.
  - `postprocess` can be used to arbitrarilly transform the synthesized (and LSF-convolved) spectrum
    before calculating the chi2.  It should take the form `postprocess(flux, data, err)` and write
    its changes in-place to the flux array.
  - `LSF_matrix`: this can be provedided along with `synthesis_wls` in place of specifying `R` if
    you have a precomputed custom LSF matrix.
  - `synthesis_wls`: see `LSF_matrix` above. This can be a Korg.Wavelengths object or any arguments
    that can be passed to its constructor, e.g. a range or vector of ranges. Wavelengths are in Å.
  - Any additional keyword arguments will be passed to [`Korg.synthesize`](@ref) when synthesizing the
    spectra for the fit.

# Returns

A NamedTuple with the following fields:

  - `best_fit_params`: the best-fit parameters
  - `best_fit_flux`: the best-fit flux, with LSF applied, resampled, and rectified.
  - `obs_wl_mask`: a bitmask for `obs_wls` which selects the wavelengths used in the fit (i.e. those
    in the `windows`)
  - `solver_result`: the result object from `Optim.jl`
  - `trace`: a vector of NamedTuples, each of which contains the parameters at each step of the
    optimization. This is empty for single parameter fits, because the underlying solver doesn't
    supply it.
  - `covariance`: a pair `(params, Σ)` where `params` is vector of parameter name (providing an
    order), and `Σ` is an estimate of the covariance matrix of the parameters.  It is the approximate
    inverse hessian of the log likelihood at the best-fit parameter calculated by the BGFS algorithm,
    and should be interpreted with caution.

!!! tip

    This function takes a long time to compile the first time it is called. Compilation performance
    is significantly better on Julia 1.10+ than previous versions, so if you are using an older
    version of Julia, you may want to upgrade.
"""
function fit_spectrum(obs_wls, obs_flux, obs_err, linelist, initial_guesses, fixed_params=(;);
                      windows=nothing, R=nothing, LSF_matrix=nothing, synthesis_wls=nothing,
                      wl_buffer=1.0, precision=1e-4, postprocess=Returns(nothing),
                      time_limit=10_000, adjust_continuum=false, synthesis_kwargs...)
    if length(obs_wls) != length(obs_flux) || length(obs_wls) != length(obs_err)
        throw(ArgumentError("obs_wls, obs_flux, and obs_err must all have the same length."))
    end
    # wavelengths, windows and LSF
    synthesis_wls, obs_wl_mask, LSF_matrix = _setup_wavelengths_and_LSF(obs_wls, synthesis_wls,
                                                                        LSF_matrix, R, windows,
                                                                        wl_buffer)

    initial_guesses, fixed_params = validate_params(initial_guesses, fixed_params)
    ps = collect(pairs(scale(initial_guesses)))
    params_to_fit = first.(ps)
    p0 = last.(ps) # the initial guess as a vector of scaled values

    @assert length(initial_guesses)>0 "Must specify at least one parameter to fit."

    chi2 = let obs_flux = obs_flux[obs_wl_mask], obs_err = obs_err[obs_wl_mask],
        obs_wls = obs_wls[obs_wl_mask], synthesis_wls = synthesis_wls, LSF_matrix = LSF_matrix,
        linelist = linelist, params_to_fit = params_to_fit, fixed_params = fixed_params

        function chi2(scaled_p)
            # this extremely weak prior helps to regularize the optimization
            negative_log_scaled_prior = sum(@. scaled_p^2 / 100^2)
            guess = unscale(Dict(params_to_fit .=> scaled_p))
            params = merge(guess, fixed_params)
            flux = postprocessed_synthetic_spectrum(synthesis_wls, linelist, LSF_matrix, params,
                                                    synthesis_kwargs, obs_wls, windows, obs_flux,
                                                    obs_err, postprocess, adjust_continuum)
            sum(((flux .- obs_flux) ./ obs_err) .^ 2) + negative_log_scaled_prior
        end
    end

    # call optimization library
    res = optimize(chi2, p0, BFGS(; linesearch=LineSearches.BackTracking(; maxstep=1.0)),
                   Optim.Options(; x_tol=precision, time_limit=time_limit, store_trace=true,
                                 extended_trace=true); autodiff=:forward)

    best_fit_params = unscale(Dict(params_to_fit .=> res.minimizer))

    best_fit_flux = try
        full_solution = merge(best_fit_params, fixed_params)
        postprocessed_synthetic_spectrum(synthesis_wls, linelist, LSF_matrix, full_solution,
                                         synthesis_kwargs, obs_wls[obs_wl_mask], windows,
                                         obs_flux[obs_wl_mask], obs_err[obs_wl_mask],
                                         postprocess, adjust_continuum)
    catch e
        println(stderr, "Exception while synthesizing best-fit spectrum")
        rethrow(e)
    end

    trace = map(res.trace) do t
        unscaled_params = unscale(Dict(params_to_fit .=> t.metadata["x"]))
        unscaled_params["chi2"] = t.value
        unscaled_params
    end

    invH = let
        # derivate relating the scaled parameters to the unscaled parameters
        # (used to convert the approximate hessian to a covariance matrix in the unscaled params)
        dp_dscaledp = map(res.minimizer, params_to_fit) do scaled_param, param_name
            ForwardDiff.derivative(scaled_param) do scaled_param
                unscale(Dict(param_name => scaled_param))[param_name]
            end
        end
        # the fact that the scaling is a diagonal operation means that we can do this as an
        # element-wise product.  If we think of ds/dp as a (diagonal) matrix, this is equivalent to
        # (ds/dp)^T * invH * (ds/dp)
        invH_scaled = res.trace[end].metadata["~inv(H)"]
        invH_scaled .* dp_dscaledp .* dp_dscaledp'
    end

    (best_fit_params=best_fit_params, best_fit_flux=best_fit_flux, obs_wl_mask=obs_wl_mask,
     solver_result=res, trace=trace, covariance=(params_to_fit, invH))
end

"""
    linear_continuum_adjustment!(obs_wls, windows, model_flux, obs_flux, obs_err)

Adjust the model flux to match the observed flux by fitting a line (as a function of wavelength) to
the residuals, and dividing it out. This can compensate for poorly done continuum normalization.

Note, obs_wls must be masked.
"""
function linear_continuum_adjustment!(obs_wls, windows, model_flux, obs_flux, obs_err)
    @assert length(obs_wls) == length(model_flux) == length(obs_flux) == length(obs_err)
    # this calculates ranges which are the same for each iteration.  Ideally, we would do this
    # ahead of time in _setup_wavelengths_and_LSF.

    if isnothing(windows)
        windows = [(first(obs_wls), last(obs_wls))]
    end

    for (λstart, λstop) in windows
        lb = searchsortedfirst(obs_wls, λstart)
        ub = searchsortedlast(obs_wls, λstop)
        ivar = 1 ./ obs_err[lb:ub] .^ 2

        X = [model_flux[lb:ub] model_flux[lb:ub] .* obs_wls[lb:ub]]
        β = (X' * (ivar .* X)) \ (X' * (ivar .* obs_flux[lb:ub]))

        view(model_flux, lb:ub) .*= β[1] .+ β[2] .* obs_wls[lb:ub]
    end
    return
end

# called by fit_spectrum
function _setup_wavelengths_and_LSF(obs_wls, synthesis_wls, LSF_matrix, R, windows, wl_buffer)
    if (!isnothing(LSF_matrix) || !isnothing(synthesis_wls))
        if !isnothing(R)
            throw(ArgumentError("LSF_matrix and synthesis_wls cannot be specified if R is provided."))
        end
        if !isnothing(windows)
            throw(ArgumentError("LSF_matrix and synthesis_wls can't be specified if windows is also specified. Passing an LSF_matrix directly is no longer needed for performance."))
        end
        if isnothing(LSF_matrix) || isnothing(synthesis_wls)
            throw(ArgumentError("LSF_matrix and synthesis_wls must both be defined if the other is."))
        end

        # do this after verifying that it's not nothing, but before checking the lengths
        synthesis_wls = Korg.Wavelengths(synthesis_wls)

        if length(obs_wls) != size(LSF_matrix, 1)
            throw(ArgumentError("the first dimension of LSF_matrix ($(size(LSF_matrix, 1))) must be the length of obs_wls ($(length(obs_wls)))."))
        end
        if (length(synthesis_wls) != size(LSF_matrix, 2))
            throw(ArgumentError("the second dimension of LSF_matrix $(size(LSF_matrix, 2)) must be the length of synthesis_wls ($(length(synthesis_wls)))"))
        end

        synthesis_wls, ones(Bool, length(obs_wls)), LSF_matrix
    else
        if isnothing(R)
            throw(ArgumentError("The resolution, R, must be specified with a keyword argument to " *
                                "Korg.Fit.fit_spectrum (unless LSF_matrix and synthesis_wls are provided.)"))
        end

        # wavelengths and windows setup
        isnothing(windows) && (windows = [(first(obs_wls), last(obs_wls))])

        windows, _ = Korg.merge_bounds(windows, 2wl_buffer)
        synthesis_wls = Korg.Wavelengths([(w[1] - wl_buffer, w[2] + wl_buffer) for w in windows])
        obs_wl_mask = zeros(Bool, length(obs_wls))
        for (λstart, λstop) in windows
            lb = searchsortedfirst(obs_wls, λstart)
            ub = searchsortedlast(obs_wls, λstop)
            obs_wl_mask[lb:ub] .= true
        end
        LSF_matrix = Korg.compute_LSF_matrix(synthesis_wls, obs_wls[obs_wl_mask], R)

        synthesis_wls, obs_wl_mask, LSF_matrix
    end
end

"""
    calculate_EWs(atm, linelist, A_X; kwargs...)

Calculate the equivalent widths of the lines in `linelist` in the spectrum synthesized from `atm`
with abundances `A_X`.

# Arguments:

  - `atm`: the model atmosphere (see [`Korg.read_model_atmosphere`](@ref) and
    [`Korg.interpolate_marcs`](@ref)).
  - `linelist`: A vector of [`Korg.Line`](@ref)s (see [`Korg.read_linelist`](@ref)).  The lines must
    be sorted by wavelength.
  - `A_X`: a vector containing the A(X) abundances (log(n_X/n_H) + 12) for elements from hydrogen to
    uranium (see [`Korg.format_A_X`](@ref)). All syntheses are done with these abundances, so if the
    resulting abundances deviate significantly from these, you may wish to iterate.

# Keyword arguments:

  - `ew_window_size` (default: 2): the farthest (in Å) to consider equivalent width contributions for
    each line.
  - `wl_step` (default: 0.01): the resolution in Å at which to synthesize the spectrum around each
    line.
  - `blend_warn_threshold` (default: 0.01): the minimum depth between two lines allowed before
    triggering a warning that they may be blended.
"""
function calculate_EWs(atm, linelist, A_X; ew_window_size::Real=2.0, wl_step=0.01,
                       blend_warn_threshold=0.01, synthesize_kwargs...)
    if !issorted(linelist; by=l -> l.wl)
        throw(ArgumentError("linelist must be sorted"))
    end

    merged_windows, lines_per_window = Korg.merge_bounds([(line.wl * 1e8 - ew_window_size,
                                                           line.wl * 1e8 + ew_window_size)
                                                          for line in linelist], 0.0)
    wl_ranges = map(merged_windows) do (wl1, wl2)
        wl1:wl_step:wl2
    end

    # hydrogen_lines should be disabled for most accurate equivalent widths.  This can be overridden
    # by passing hydrogen_lines=true as a keyword argument (included in synthesize_kwargs)
    # line_buffer=0.0 makes things a bit faster, and it causes no problems as long as ew_window_size
    # is sufficient, which is necessary anyway.
    sol = Korg.synthesize(atm, linelist, A_X, wl_ranges; line_buffer=0.0, hydrogen_lines=false,
                          synthesize_kwargs...)
    depth = 1 .- sol.flux ./ sol.cntm

    element_type = promote_type(eltype(A_X), eltype(Korg.get_temps(atm)),
                                typeof(linelist[1].log_gf)) #TODO hacky
    EWs = Array{element_type}(undef, length(linelist))
    all_boundaries = Float64[]
    for (wl_range, subspec, line_indices) in zip(wl_ranges, sol.subspectra, lines_per_window)
        absorption = depth[subspec]

        # get the wl-index of least absorption between each pair of lines
        boundary_indices = map(1:length(line_indices)-1) do i
            wl1 = linelist[line_indices[i]].wl * 1e8
            wl2 = linelist[line_indices[i+1]].wl * 1e8
            l1_ind = Int(round((wl1 - wl_range[1]) / step(wl_range))) + 1
            l2_ind = Int(round((wl2 - wl_range[1]) / step(wl_range))) + 1
            boundary_index = argmin(absorption[l1_ind:l2_ind]) + l1_ind - 1
            if absorption[boundary_index] > blend_warn_threshold
                @warn "Lines $(line_indices[i]) and $(line_indices[i+1]) ($(linelist[line_indices[i]].wl*1e8) Å and $(linelist[line_indices[i+1]].wl*1e8)) Å appear to be blended.  Between them, the absorption never drops below $(blend_warn_threshold) (minimum: $(ForwardDiff.value(absorption[boundary_index]))). You can adjust this threshold with the blend_warn_threshold keyword argument."
            end
            boundary_index
        end
        boundary_indices = [1; boundary_indices; length(subspec)]
        for b in boundary_indices
            push!(all_boundaries, wl_range[b])
        end

        for i in 1:length(line_indices)
            r = boundary_indices[i]:boundary_indices[i+1]
            EWs[line_indices[i]] = trapz(wl_range[r], absorption[r]) * 1e3 # convert to mÅ
        end
    end
    EWs
end

"""
    ews_to_abundances(atm, linelist, A_X, measured_EWs; kwargs... )

!!! danger

    Functions using equivalent widths have been dissabled while we
    [fix some bugs](https://github.com/ajwheeler/Korg.jl/pull/331).

Compute per-line abundances on the linear part of the curve of growth given a model atmosphere and a
list of lines with equivalent widths.

TODO assumptions

# Arguments:

  - `atm`: the model atmosphere (see [`Korg.read_model_atmosphere`](@ref) and
    [`Korg.interpolate_marcs`](@ref)).
  - `linelist`: A vector of [`Korg.Line`](@ref)s (see [`Korg.read_linelist`](@ref)).  The lines must
    be sorted by wavelength.
  - `A_X`: a vector containing the A(X) abundances (log(n_X/n_H) + 12) for elements from hydrogen to
    uranium (see [`Korg.format_A_X`](@ref)). All syntheses are done with these abundances, so if the
    resulting abundances deviate significantly from these, you may wish to iterate.
  - `measured_EWs`: a vector of equivalent widths (in mÅ)

# Returns

A vector of abundances (`A(X) = log10(n_X/n_H) + 12` format) for each line in `linelist`.

# Optional arguments:

  - `wl_step` (default: 0.01) is the resolution in Å at which to synthesize the spectrum around each
    line.
  - `ew_window_size` (default: 2): the farthest (in Å) to consider equivalent width contributions for
    each line.  It's very important that this is large enough to include each line entirely.
  - `blend_warn_threshold` (default: 0.01) is the minimum absorption between two lines allowed before
    triggering a warning that they may be blended.
    All other keyword arguments are passed to [`Korg.synthesize`](@ref) when synthesizing each line.
  - `finite_difference_delta_A` (default: 0.01): the step size in A(X) to use for the finite
    difference calculation of the curve of growth slope.
"""
# TODO

"""
    assume that all lines are on the linear part of the curve of growth?
"""
function ews_to_abundances_approx(atm, linelist, A_X, measured_EWs; ew_window_size::Real=2.0,
                                  wl_step=0.01,
                                  blend_warn_threshold=0.01, finite_difference_delta_A=0.01,
                                  synthesize_kwargs...)
    if length(linelist) != length(measured_EWs)
        throw(ArgumentError("length of linelist does not match length of ews ($(length(linelist)) != $(length(measured_EWs)))"))
    end

    if any(l -> Korg.ismolecule(l.species), linelist)
        throw(ArgumentError("linelist contains molecular species"))
    end

    # Check that the user is supplying EWs in mA
    if 1 > maximum(measured_EWs)
        @warn "Maximum EW given is less than 1 mA. Check that you're giving EWs in mÅ (*not* Å)."
    end

    A0 = [A_X[Korg.get_atoms(l.species)[1]] for l in linelist]

    EWs = calculate_EWs(atm, linelist, A_X; ew_window_size=ew_window_size, wl_step=wl_step,
                        blend_warn_threshold=blend_warn_threshold, synthesize_kwargs...)
    # ∂A/∂REW = 1 by assumption
    @. A0 + (log10(measured_EWs) - log10.(EWs))
end

"""
TODO

TODO rename
"""
function ews_to_abundances(atm, linelist, A_X, measured_EWs; ew_window_size::Real=2.0,
                           wl_step=0.01, callback=Returns(nothing), abundance_tol=1e-5,
                           max_iter=30, blend_warn_threshold=0.01, synthesize_kwargs...)
    # TODO factor out these checks if multiple methods remain
    if length(linelist) != length(measured_EWs)
        throw(ArgumentError("length of linelist does not match length of ews ($(length(linelist)) != $(length(measured_EWs)))"))
    end

    if any(l -> Korg.ismolecule(l.species), linelist)
        throw(ArgumentError("linelist contains molecular species"))
    end

    # Check that the user is supplying EWs in mA
    if 1 > maximum(measured_EWs)
        @warn "Maximum EW given is less than 1 mA. Check that you're giving EWs in mÅ (*not* Å)."
    end

    # do a single synthesis to get the chemical equilibrium once
    sol = synthesize(atm, [], A_X, 5000, 5000)

    A_X = copy(A_X)
    # fiducial abundance for each line is taked from the A_X vector
    A0 = [A_X[Korg.get_atom(l.species)] for l in linelist]
    # difference between the measured and fiducial abundance for each line
    ΔA = zeros(eltype(sol.flux), length(linelist))
    ΔA_prev = copy(ΔA) .+ 0.01
    fitmask = ones(Bool, length(linelist))

    perturbed_linelist = [Korg.Line(l; log_gf=l.log_gf .+ Δ)
                          for (l, Δ) in zip(linelist, ΔA_prev - ΔA)]

    EWs_prev = calculate_EWs(atm, perturbed_linelist, A_X; ew_window_size=ew_window_size,
                             wl_step=wl_step, blend_warn_threshold=blend_warn_threshold,
                             use_chemical_equilibrium_from=sol, synthesize_kwargs...)
    EWs = copy(EWs_prev) # just to allocate

    iter = 0
    while sum(fitmask) > 0 && iter < max_iter # while there are still lines to fit
        iter += 1
        perturbed_linelist = [Korg.Line(l; log_gf=l.log_gf .+ Δ)
                              for (l, Δ) in zip(linelist[fitmask], ΔA[fitmask])]
        EWs[fitmask] .= calculate_EWs(atm, perturbed_linelist, A_X; ew_window_size=ew_window_size,
                                      wl_step=wl_step, blend_warn_threshold=blend_warn_threshold,
                                      use_chemical_equilibrium_from=sol, synthesize_kwargs...)

        ∂A_∂REW = @. (ΔA[fitmask] - ΔA_prev[fitmask]) / log10(EWs[fitmask] / EWs_prev[fitmask])
        δA = @. ∂A_∂REW * log10(measured_EWs[fitmask] / EWs[fitmask])

        ΔA_prev[fitmask] .= ΔA[fitmask]
        ΔA[fitmask] .+= δA
        callback(A0 .+ ΔA)

        # stop fitting lines that have converged
        fitmask[fitmask] .&= abs.(δA) .> abundance_tol

        # TODO verbose?
        println("iter $iter ($(sum(fitmask)) lines unconverged)")
    end

    ΔA[fitmask] .= NaN # NaN out anything unconverged

    A0 .+ ΔA
end

"""
    _validate_stellar_parameters(linelist, measured_EWs, measured_EW_err, params0, parameter_ranges)

Validate the input parameters for stellar parameter determination.
"""
function _validate_stellar_parameters(linelist, measured_EWs, measured_EW_err, params0,
                                      parameter_ranges)
    if length(linelist) != length(measured_EWs) || length(linelist) != length(measured_EW_err)
        throw(ArgumentError("length of linelist does not match length of ews ($(length(linelist)) != $(length(measured_EWs)))"))
    end

    formulas = [line.species.formula for line in linelist]
    if any(Ref(formulas[1]) .!= formulas)
        throw(ArgumentError("All lines must be from the same element."))
    end

    if Korg.ismolecule(linelist[1].species)
        throw(ArgumentError("Cannot do stellar parameter determination with molecular lines."))
    end

    neutrals = [l.species.charge == 0 for l in linelist]
    if (sum(neutrals) < 3) || (sum(.!neutrals) < 1)
        throw(ArgumentError("Must have at least 3 neutral lines and 1 ion line."))
    end

    if params0[3] == 0.0  # vmic0
        throw(ArgumentError("Starting guess for vmic (vmic0) must be nonzero."))
    end

    if any(p[1] >= p[2] for p in parameter_ranges)
        throw(ArgumentError("The lower bound of each parameter must be less than the upper bound."))
    end

    if parameter_ranges[3][1] <= 0.0
        throw(ArgumentError("The lower bound of vmic must be greater than zero. (vmic must be nonzero in order to avoid null derivatives. Very small values are fine.)"))
    end

    # the widest parameter ranges allowed for model atmosphere interp
    atm_lb = first.(Korg._sdss_marcs_atmospheres[1][1:3])
    atm_lb[3] = Korg._low_Z_marcs_atmospheres[1][3][1]
    atm_ub = last.(Korg._sdss_marcs_atmospheres[1][1:3])
    if any(first.(parameter_ranges[[1, 2, 4]]) .< atm_lb) ||
       any(last.(parameter_ranges[[1, 2, 4]]) .> atm_ub)
        throw(ArgumentError("The parameter ranges must be within the range of the MARCS grid"))
    end

    params = clamp(params0, first.(parameter_ranges), last.(parameter_ranges))
    for (p, p0, n) in zip(params, params0, ["Teff", "logg", "vmic", "metallicity"])
        if p != p0
            @warn "Initial guess for $n ($p0) has been clamped to $p, to be within the allowed range."
        end
    end

    return params
end

"""
    ews_to_stellar_parameters(linelist, measured_EWs, [measured_EW_err]; kwargs...)

!!! danger

    Functions using equivalent widths have been dissabled while we
    [fix some bugs](https://github.com/ajwheeler/Korg.jl/pull/331).

Find stellar parameters from equivalent widths the "old fashioned" way.  This function finds the
values of ``T_\\mathrm{eff}``, ``\\log g``, ``v_{mic}``, and [m/H] which satisfy the following conditions
(using a Newton-Raphson solver):

  - The slope of the abundances of neutral lines with respect to lower excitation potential is zero.
  - The difference between the mean abundances of neutral and ionized lines is zero.
  - The slope of the abundances of neutral lines with respect to reduced equivalent width is zero.
  - The difference between the mean abundances of all lines and the model-atmosphere input [m/H] is zero.
    Here the "slope" refers to the slope of a linear fit to the abundances of the lines in question.

# Arguments:

  - `linelist`: A vector of [`Korg.Line`](@ref) objects (see [`Korg.read_linelist`](@ref)).  The lines
    must be sorted by wavelength.
  - `measured_EWs`: a vector of equivalent widths (in mÅ).
  - `measured_EW_err` (optional): the uncertainty in `measured_EWs`.  If not specified, all lines are
    assumed to have the same uncertainty. These uncertainties are used when evaluating the equations
    above, and are propagated to provide uncertainties in the resulting parameters.

# Returns:

A tuple containing:

  - the best-fit parameters: `[Teff, logg, vmic, [m/H]]` as a vector
  - the statistical uncertainties in the parameters, propagated from the uncertainties in the
    equivalent widths. This is zero when EW uncertainties are not specified.
  - the systematic uncertainties in the parameters, estimated from the scatter the abundances computed
    from each line not accounted for by the EW uncertainties.

# Keyword arguments:

  - `Teff0` (default: 5000.0) is the starting guess for Teff

  - `logg0` (default: 3.5) is the starting guess for logg
  - `vmic0` (default: 1.0) is the starting guess for vmic. Note that this must be nonzero in order to
    avoid null derivatives. Very small values are fine.
  - `m_H0` (default: 0.0) is the starting guess for [m/H]
  - `tolerances` (default: `[1e-3, 1e-3, 1e-4, 1e-3]`) is the tolerance for the residuals each equation
    listed above. The solver stops when all residuals are less than the corresponding tolerance.
  - `max_step_sizes` (default: `[1000.0, 1.0, 0.3, 0.5]`) is the maximum step size to take in each
    parameter direction.  This is used to prevent the solver from taking too large of a step and
    missing the solution.  Be particularly cautious with the vmic (third) parameter, as the
    unadjusted step size is often too large.
  - `parameter_ranges` (default: `[(2800.0, 8000.0), (-0.5, 5.5), (1e-3, 10.0), (-2.5, 1.0)]`) is the
    allowed range for each parameter. This is used to prevent the solver from wandering into
    unphysical parameter space, or outside the range of the MARCS grid supported by
    [`Korg.interpolate_marcs`](@ref). The default ranges ``T_\\mathrm{eff}``, ``\\log g``,  and [m/H]
    are the widest supported by the MARCS grid. Note that vmic must be nonzero in order to avoid null
    derivatives.
  - `callback`:  is a function which is called at each step of the optimization.
    It is passed three arguments:

      + the current values of the parameters
      + the residuals of each equation being solved
      + the abundances of each line computed with the current parameters.
        You can pass a callback function, to e.g. make a plot of the residuals at each step.
  - `max_iterations` (default: 30) is the maximum number of iterations to allow before stopping the
    optimization.
"""
function ews_to_stellar_parameters(linelist, measured_EWs,
                                   measured_EW_err=ones(length(measured_EWs));
                                   Teff0=5000.0, logg0=3.5, vmic0=1.0, m_H0=0.0,
                                   tolerances=[1e-3, 1e-3, 1e-4, 1e-3],
                                   max_step_sizes=[1000.0, 1.0, 0.3, 0.5],
                                   # TODO expand to include low metallicity atmospheres
                                   parameter_ranges=[(2800.0, 8000.0),
                                       (-0.5, 5.5),
                                       (1e-3, 10.0),
                                       (-2.5, 1.0)],
                                   fix_params=[false, false, false, false],
                                   callback=Returns(nothing), max_iterations=30, passed_kwargs...)
    if :vmic in keys(passed_kwargs)
        throw(ArgumentError("vmic must not be specified, because it is a parameter fit by ews_to_stellar_parameters.  Did you mean to specify vmic0, the starting value? See the documentation for ews_to_stellar_parameters if you would like to fix microturbulence to a given value."))
    end

    params0 = [Teff0, logg0, vmic0, m_H0]
    params = _validate_stellar_parameters(linelist, measured_EWs, measured_EW_err, params0,
                                          parameter_ranges)

    # set up closure to compute residuals
    get_residuals = (p) -> _stellar_param_equation_residuals(false, p, linelist, measured_EWs,
                                                             measured_EW_err, fix_params, callback,
                                                             passed_kwargs)
    iterations = 0
    J_result = DiffResults.JacobianResult(params)
    while true
        _ews_to_stellar_parameters_iteration!(J_result, get_residuals, params, fix_params,
                                              tolerances, max_step_sizes, parameter_ranges) && break
        iterations += 1
        if iterations > max_iterations
            @warn "Failed to converge after $max_iterations iterations.  Returning the current guess."
            return params, fill(NaN, 4), fill(NaN, 4)
        end
    end
    iterations = 0
    get_residuals = (p) -> _stellar_param_equation_residuals(true, p, linelist, measured_EWs,
                                                             measured_EW_err, fix_params, callback,
                                                             passed_kwargs)
    while true
        _ews_to_stellar_parameters_iteration!(J_result, get_residuals, params, fix_params,
                                              tolerances, max_step_sizes, parameter_ranges) && break
        iterations += 1
        if iterations > max_iterations
            @warn "Failed to converge after $max_iterations iterations.  Returning the current guess."
            return params, fill(NaN, 4), fill(NaN, 4)
        end
    end

    # compute uncertainties
    stat_σ_r, sys_σ_r = _stellar_param_residual_uncertainties(params, linelist, measured_EWs,
                                                              measured_EW_err, passed_kwargs)

    J = DiffResults.jacobian(J_result)[.!fix_params, .!fix_params]
    stat_σ = zeros(4)
    stat_σ[.!fix_params] .= abs.(J \ stat_σ_r[.!fix_params])
    sys_σ = zeros(4)
    sys_σ[.!fix_params] .= abs.(J \ sys_σ_r[.!fix_params])

    params, stat_σ, sys_σ
end

function _ews_to_stellar_parameters_iteration!(J_result, get_residuals, params,
                                               fix_params, tolerances, max_step_sizes,
                                               parameter_ranges)
    J_result = ForwardDiff.jacobian!(J_result, get_residuals, params)
    J = DiffResults.jacobian(J_result)
    residuals = DiffResults.value(J_result)
    if all((abs.(residuals).<tolerances)[.!fix_params]) # stopping condition
        return true
    end

    # calculate step, and update params
    step = zeros(length(params))
    step[.!fix_params] = -J[.!fix_params, .!fix_params] \ residuals[.!fix_params]
    params .+= clamp.(step, -max_step_sizes, max_step_sizes)
    params .= clamp.(params, first.(parameter_ranges), last.(parameter_ranges))

    return false
end

# called by ews_to_stellar_parameters
function _stellar_param_equation_residuals(exact_calculation, params, linelist, EW, EW_err,
                                           fix_params, callback,
                                           passed_kwargs)
    A, A_inv_var, neutrals, REWs, Z = _stellar_param_equations_precalculation(exact_calculation,
                                                                              params, linelist, EW,
                                                                              EW_err, passed_kwargs)

    finitemask = isfinite.(A)
    neutrals = neutrals .& finitemask

    teff_residual = get_slope([line.E_lower for line in linelist[neutrals]],
                              A[neutrals], A_inv_var[neutrals])
    logg_residual = (weighted_mean(A[neutrals], A_inv_var[neutrals]) -
                     weighted_mean(A[.!neutrals.&finitemask], A_inv_var[.!neutrals.&finitemask]))

    # TODO make sure it works if a line fails to converge
    #@show length(REWs)
    #@show length(finitemask)
    #@show length(neutrals)
    #@show length(finitemask[neutrals])

    vmic_residual = get_slope(REWs[finitemask[neutrals]], A[neutrals], A_inv_var[neutrals])
    feh_residual = weighted_mean(A[finitemask], A_inv_var[finitemask]) -
                   (params[4] + Korg.grevesse_2007_solar_abundances[Z])
    residuals = [teff_residual, logg_residual, vmic_residual, feh_residual]
    residuals .*= .!fix_params # zero out residuals for fixed parameters

    callback(ForwardDiff.value.(params), ForwardDiff.value.(residuals), ForwardDiff.value.(A))
    residuals
end

# called by _stellar_param_equation_residuals
# returns (statistical_uncertainty, systematic_uncertainty)
function _stellar_param_residual_uncertainties(params, linelist, EW, EW_err, passed_kwargs)
    A, A_inv_var, neutrals, REWs, _ = _stellar_param_equations_precalculation(true, params,
                                                                              linelist, EW,
                                                                              EW_err, passed_kwargs)

    # estimated total (including systematic) err in the abundances of each line
    total_err = std(A)
    total_ivar = ones(length(A)) * total_err^-2

    stat_sigma, total_sigma = map([A_inv_var, total_ivar]) do ivar
        sigma_mean = 1 ./ sqrt(sum(ivar))
        teff_residual_sigma = get_slope_uncertainty([line.E_lower for line in linelist[neutrals]],
                                                    ivar[neutrals])
        vmic_residual_sigma = get_slope_uncertainty(REWs, ivar[neutrals])
        [teff_residual_sigma, sigma_mean, vmic_residual_sigma, sigma_mean]
    end

    if EW_err == ones(length(EW))
        # in the case that EW uncertainties were specified, the residuals were calculated
        # assuming errs = 1, but we want to ignore them here, and call all error "systematic"
        [0.0, 0.0, 0.0, 0.0], total_sigma
    else
        stat_sigma, sqrt.(max.(total_sigma .^ 2 .- stat_sigma .^ 2, 0))
    end
end

function _stellar_param_equations_precalculation(exact_calculation, params, linelist, EW, EW_err,
                                                 passed_kwargs)
    teff, logg, vmic, m_H = params
    A_X = Korg.format_A_X(m_H)
    atm = Korg.interpolate_marcs(teff, logg, A_X; perturb_at_grid_values=true,
                                 clamp_abundances=true)

    if exact_calculation
        A = ews_to_abundances(atm, linelist, A_X, EW; vmic=vmic, blend_warn_threshold=Inf,
                              passed_kwargs...)
    else
        A = ews_to_abundances_approx(atm, linelist, A_X, EW; vmic=vmic, blend_warn_threshold=Inf,
                                     passed_kwargs...)
    end
    # convert error in EW to inverse variance in A (assuming linear part of C.O.G.)
    A_inv_var = (EW .* A ./ EW_err) .^ 2

    neutrals = [l.species.charge == 0 for l in linelist]
    REWs = log10.(EW[neutrals] ./ [line.wl for line in linelist[neutrals]])
    # this is guaranteed not to be a mol (checked by ews_to_stellar_parameters).
    Z = Korg.get_atoms(linelist[1].species)[1]

    A, A_inv_var, neutrals, REWs, Z
end

weighted_mean(x, inv_var) = sum(x .* inv_var) / sum(inv_var)

# called by _stellar_param_equation_residuals
function get_slope(xs, ys, inv_var)
    Δx = xs .- mean(xs)
    Δy = ys .- mean(ys)
    sum(Δx .* Δy .* inv_var) ./ sum(Δx .^ 2 .* inv_var)
end

function get_slope_uncertainty(xs, ivar::AbstractVector) # guard against scalar ivar
    sqrt(sum(ivar) / (sum(ones(length(xs)) .* ivar) * sum(ivar .* xs .^ 2) - sum(ivar .* xs)^2))
end

end # module
