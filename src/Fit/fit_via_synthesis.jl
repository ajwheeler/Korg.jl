using LsqFit
using LinearAlgebra: cond
using Interpolations: linear_interpolation, Line

# Bounds for each paramter the fit_spectrum can fit. These are passed to the Levenberg-Marquardt 
# optimizer, and also define the linear scaling that the optimization works in (see scale_param).
const param_bounds = Dict("epsilon" => (0.0, 1.0),
                          "cntm_offset" => (-0.5, 0.5),
                          "cntm_slope" => (-0.1, 0.1),
                          # These limits come from what is supported by Korg.interpolate_marcs.
                          "Teff" => (2800.0, 8000.0),
                          "logg" => (-0.5, 5.5),
                          "M_H" => (-5.0, 1.0),
                          # this allows all the atmospheres supported by the grid, but also many that are not.
                          # alpha will be clamped to the nearest supported value.
                          "alpha_H" => (-3.5, 2.0),
                          # vmic and vsini must be non-negative. The lower bound on vmic is a small
                          # positive number to avoid null derivatives.
                          "vmic" => (1e-3, 250.0),
                          "vsini" => (0.0, 250.0),
                          map(Korg.atomic_symbols) do el
                              el => (-10.0, +4.0)
                          end...)

# The optimizer works in scaled coordinates, in which each parameter varies between 0 and 1. This is 
# to make the parameters comparably sized, so that convergence conditions in terms of x_tol are 
# sensible.
function scale_param(p, name)
    (lower, upper) = param_bounds[name]
    (p - lower) / (upper - lower)
end
function unscale_param(p, name)
    (lower, upper) = param_bounds[name]
    p * (upper - lower) + lower
end

"""
Synthesize a spectrum, returning the flux, with LSF applied, resampled, and rectified.  This is
an internal function oned by fitting routines.
See [`Korg.synthesize`](@ref) to synthesize spectra as a Korg user.
"""
function synthetic_spectrum(synthesis_wls, linelist, LSF_matrix, params, synthesis_kwargs)
    specified_abundances = Dict([p for p in pairs(params) if p.first in Korg.atomic_symbols])
    alpha_H = "alpha_H" in keys(params) ? params["alpha_H"] : params["M_H"]
    A_X::Vector{valtype(params)} = Korg.format_A_X(params["M_H"], alpha_H, specified_abundances;
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
        if (e isa Korg.ChemicalEquilibriumError) || (e isa Korg.LazyMultilinearInterpError) ||
           (e isa Korg.AtmosphereInterpolationError)
            # These errors happen at a few parameter values that synthesis can't handle.
            # LazyMultilinearInterpError: happens when logg is out of bounds for the low-Z 
            # atmosphere grid (can't be expressed with box constraints.) The other two are rarer 
            # failure modes.
            # Return Inf as the residuals at all pixels, to ensure the step is rejected.
            return fill(oftype(float(first(obs_err)), Inf), length(obs_flux))
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
                         default_params=Dict("M_H" => 0.0, "vsini" => 0.0, "vmic" => 1.0,
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
        if param == "m_H"
            throw(ArgumentError("m_H is no longer a supported keyword argument of fit_spectrum (starting in Korg 1.0). Use M_H instead."))
        end
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

    for (param, value) in [pairs(initial_guesses)...; pairs(fixed_params)...]
        (lower, upper) = param_bounds[param]
        if !(lower <= value <= upper)
            throw(ArgumentError("$param = $value is outside the range supported by " *
                                "Korg.Fit.fit_spectrum, $lower to $upper."))
        end
    end

    initial_guesses, fixed_params
end

# make it possible to use dicts instead of NamedTuples for the python people
function validate_params(initial_guesses::AbstractDict, fixed_params::NamedTuple; kwargs...)
    validate_params(initial_guesses, _namedtuple_to_dict(fixed_params); kwargs...)
end
function validate_params(initial_guesses::NamedTuple, fixed_params; kwargs...)
    validate_params(_namedtuple_to_dict(initial_guesses), fixed_params; kwargs...)
end

function _namedtuple_to_dict(nt::NamedTuple)
    Dict{String,Float64}([string(p[1]) => Float64(p[2]) for p in pairs(nt)])
end

# throw informative errors if the observed spectrum has issues
function _validate_observed_spectrum(obs_wls, obs_flux, obs_err, obs_wl_mask)
    if length(obs_wls) != length(obs_flux) || length(obs_wls) != length(obs_err)
        throw(ArgumentError("When using Korg.Fit.fit_spectrum, obs_wls, obs_flux, and obs_err must all have the same length."))
    end

    for arr in [obs_wls, obs_flux, obs_err]
        if any(.!isfinite, arr[obs_wl_mask])
            throw(ArgumentError("When using Korg.Fit.fit_spectrum, obs_wls, obs_flux, and obs_err must not contain NaN or Inf."))
        end
    end

    if any(iszero, obs_err[obs_wl_mask])
        throw(ArgumentError("When using Korg.Fit.fit_spectrum, obs_err must not contain zeros."))
    end
end

"""
    fit_spectrum(obs_wls, obs_flux, obs_err, linelist, initial_guesses, fixed_params; kwargs...)

Find the parameters and abundances that best match a rectified (continuum-normalized) observed
spectrum.

# Arguments:

  - `obs_wls`: the wavelengths of the observed spectrum, in any format accepted by synthesize
    (see [Wavelengths](https://ajwheeler.github.io/Korg.jl/stable/Wavelengths/))
  - `obs_flux`: the observed flux
  - `obs_err`: the uncertainty in the observed flux
  - `linelist`: the linelist to use for synthesis
  - `initial_guesses`: a NamedTuple containing initial guesses for the parameters to fit.  See
    "Specifying parameters" below.
  - `fixed_params`: a NamedTuple containing parameters to hold fixed during fitting (default: empty).
    See "Specifying parameters" below.

# Specifying parameters

`initial_guesses` and `fixed_params` can also be specified as Dicts instead of NamedTuples, which is
more convenient when calling Korg from python.

Parameters are specified as named tuples or dictionaries. Named tuples look like this:
`(Teff=5000, logg=4.5, M_H=0.0)`.  Single-element named tuples require a semicolon: `(; Teff=5000)`.

### Required parameters

`Teff` and `logg` *must* be specified in either `initial_guesses` or `fixed_params`.

### Optional Parameters

These can be specified in either `initial_guesses` or `fixed_params`, but if they are not default
values are used.

  - `M_H`: the metallicity of the star, in dex. Default: `0.0`
  - `alpha_H`: the alpha enhancement of the star, in dex. Default: `M_H`.  Note that, because of the
    parameter range supported by [`Korg.interpolate_marcs`](@ref), only values within ±1 of `M_H`
    are supported.
  - `vmic`: the microturbulence velocity, in km/s. Default: `1.0`
  - `vsini`: the projected rotational velocity of the star, in km/s. Default: `0.0`.
    See [`Korg.apply_rotation`](@ref) for details.
  - `epsilon`: the linear limb-darkening coefficient. Default: `0.6`. Used for applying rotational
    broadening only.  See [`Korg.apply_rotation`](@ref) for details.
  - Individual elements, e.g. `Na`, specify the solar-relative ([X/H]) abundance of that element.

# Keyword arguments

  - `R`, the resolution of the observed spectrum. This is required, unless you specify `LSF_matrix`
    directly.  It can be specified as a function of wavelength, in which case it will be evaluated
    at the observed wavelengths.
  - `windows` is a vector of wavelength pairs, specifying the range(s) to jointly fit. If `windows`
    is not specified, the entire spectrum is used. Overlapping windows are automatically merged.
  - `adjust_continuum` (default: `false`) if true, adjust the continuum with the best-fit linear
    correction within each window, minimizing the chi-squared between data and model at every step
    of the optimization. Note that this will result in underestimated parameter uncertainty.
  - `wl_buffer` is the number of Å to add to each side of the synthesis range for each window.
  - `time_limit` is the maximum number of seconds to spend in the optimizer. (default: `10_000`).
    The optimizer will only checks against the time limit after each step, so the actual wall time
    may exceed this limit.
  - `precision` specifies the tolerance for the solver to accept a solution (passed to the
    Levenberg-Marquardt optimizer as `x_tol`). The solver operates on scaled parameters, in which
    each parameter's allowed range maps onto `[0, 1]`, so `precision` is roughly a fraction of each
    parameter's full range. The default value, `1e-4`, corresponds to a worst-case tolerance of about
    0.5 K in `Teff`, 0.0006 in `logg`, 0.0006 in `M_H`, and 0.001 in detailed abundances.
  - `postprocess` can be used to arbitrarilly transform the synthesized, LSF-convolved spectrum
    before calculating the chi2.  It should take the form `postprocess(flux, data, err)` and write
    its changes in-place to the flux array.
  - `condition_number_warning_threshold` (default: `1e3`): if the condition
    number of the approximation of the covariance at the solution (the `condition_number` field of
    the returned object) exceeds this value, a warning is emitted that the fit is
    poorly-constrained.
  - `LSF_matrix`: this can be provedided along with `synthesis_wls` in place of specifying `R` if
    you have a precomputed custom LSF matrix.
  - `synthesis_wls`: see `LSF_matrix` above. This can be a Korg.Wavelengths object or any arguments
    that can be passed to its constructor, e.g. a range or vector of ranges. See
    [Wavelengths](https://ajwheeler.github.io/Korg.jl/stable/Wavelengths/) for all the ways
    wavelengths can be specified. Wavelengths are in Å.
  - Any additional keyword arguments will be passed to [`Korg.synthesize`](@ref) when synthesizing the
    spectra for the fit.

# Returns

An object with the following fields:

  - `best_fit_params`: the best-fit parameters
  - `best_fit_flux`: the best-fit flux, with LSF applied, resampled, and rectified.
  - `obs_wl_mask`: a bitmask for `obs_wls` which selects the wavelengths used in the fit (i.e. those
    in the `windows`)
  - `solver_result`: the result object from `LsqFit.jl`
  - `trace`: a vector of Dicts, one per optimization step, each containing the χ² (`"chi2"`), the
    infinity norm of the gradient (`"g_norm"`), and the Levenberg-Marquardt damping parameter
    (`"lambda"`) at that step.
  - `covariance`: a pair `(params, Σ)` where `params` is vector of parameter name (providing an
    order), and `Σ` is an estimate of the covariance matrix of the parameters.  It is the approximate
    inverse hessian of the log likelihood at the best-fit parameter calculated by the optimizer,
    and should be interpreted with caution.
  - `condition_number`: the condition number of the weighted Gauss-Newton approximation to the
    Hessian at the solution, in the optimizer's scaled coordinates. A large value indicates a
    poorly-constrained, near-degenerate fit whose best-fit parameters may be unreliable even when
    χ² is small.
"""
function fit_spectrum(obs_wls, obs_flux, obs_err, linelist, initial_guesses, fixed_params=(;);
                      windows=nothing, R=nothing, LSF_matrix=nothing, synthesis_wls=nothing,
                      wl_buffer=1.0, precision=1e-4, postprocess=Returns(nothing),
                      time_limit=10_000, adjust_continuum=false,
                      condition_number_warning_threshold=1e3, synthesis_kwargs...)
    if adjust_continuum
        @warn "Note that setting adjust_continuum=true will result in underestimated best-fit parameter uncertainty."
    end

    # wavelengths, windows and LSF
    (synthesis_wls, obs_wl_mask,
    LSF_matrix) = _setup_wavelengths_and_LSF(obs_wls, synthesis_wls, LSF_matrix, R, windows,
                                             wl_buffer)

    _validate_observed_spectrum(obs_wls, obs_flux, obs_err, obs_wl_mask)

    initial_guesses, fixed_params = validate_params(initial_guesses, fixed_params)
    ps = collect(pairs(initial_guesses))
    params_to_fit = first.(ps)
    # the initial guess in scaled coords
    p0 = [scale_param(p, name) for (name, p) in ps]

    @assert length(initial_guesses)>0 "Must specify at least one parameter to fit."

    # the model that LsqFit optimizes: given scaled parameters, return the synthesized flux.
    masked_obs_flux = obs_flux[obs_wl_mask]
    masked_obs_err = obs_err[obs_wl_mask]
    masked_obs_wls = obs_wls[obs_wl_mask]
    function model(_, scaled_p)
        guess = Dict(name => unscale_param(p, name)
                     for (name, p) in zip(params_to_fit, scaled_p))
        params = merge(guess, fixed_params)
        postprocessed_synthetic_spectrum(synthesis_wls, linelist, LSF_matrix, params,
                                         synthesis_kwargs, masked_obs_wls, windows,
                                         masked_obs_flux, masked_obs_err, postprocess,
                                         adjust_continuum)
    end

    # call optimization library
    lower_bounds = zeros(length(params_to_fit))
    upper_bounds = ones(length(params_to_fit))
    result = curve_fit(model,
                       obs_wls[obs_wl_mask],
                       obs_flux[obs_wl_mask],
                       1 ./ obs_err[obs_wl_mask] .^ 2, # weights are inverse variance
                       p0;
                       lower=lower_bounds,
                       upper=upper_bounds,
                       x_tol=precision,
                       maxTime=Float64(time_limit),
                       store_trace=true,
                       autodiff=:forwarddiff)
    best_fit_params = Dict(name => unscale_param(p, name)
                           for (name, p) in zip(params_to_fit, result.param))

    # raise warnings if the fit seems bad, and get the condition number of Σ
    condition_number = sanity_check_result(result, condition_number_warning_threshold,
                                           params_to_fit)

    # recover the best-fit flux from the solver's stored residuals rather than re-synthesizing.
    # LsqFit stores resid = √wt .* (model - obs_flux) evaluated at the converged parameters, and we
    # pass wt = 1 ./ obs_err.^2, so model = resid .* obs_err + obs_flux.
    best_fit_flux = result.resid .* obs_err[obs_wl_mask] .+ obs_flux[obs_wl_mask]

    trace = map(result.trace) do t
        Dict("chi2" => t.value, "g_norm" => t.g_norm, "lambda" => t.metadata["lambda"])
    end

    # covariance of the best-fit parameters, computed from the Jacobian of the model at the solution
    scales = map(params_to_fit) do name
        l, u = param_bounds[name]
        l - u
    end
    Σ = vcov(result) .* scales .* scales'

    (best_fit_params=best_fit_params, best_fit_flux=best_fit_flux, obs_wl_mask=obs_wl_mask,
     solver_result=result, trace=trace, covariance=(params_to_fit, Σ),
     condition_number=condition_number)
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
    if !issorted(obs_wls)
        throw(ArgumentError("When using Korg.Fit.fit_spectrum, obs_wls must be sorted in order of increasing wavelength."))
    end

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

function sanity_check_result(result, condition_number_warning_threshold, params_to_fit)
    # Diagnose whether the problem is degenerate. 
    # J' J is the inverse covariance of the scaled params.
    condition_number = cond(result.jacobian' * result.jacobian)
    if condition_number > condition_number_warning_threshold
        @warn "The fit is poorly constrained (condition number of JᵀWJ is $condition_number, " *
              "above the condition_number_warning_threshold of " *
              "$condition_number_warning_threshold). The best-fit parameters may be unreliable: " *
              "different parameter combinations likely produce nearly indistinguishable spectra. " *
              "Consider fitting fewer parameters, or using more/broader wavelength windows to break " *
              "the degeneracy. (Raise condition_number_warning_threshold to silence this.)"
    end
    # flag parameters that ended up pinned against their box bounds, which also indicates a bad fit
    pinned = [name for (name, p) in zip(params_to_fit, result.param)
              if p < 1e-4 || p > 1 - 1e-4]
    if !isempty(pinned)
        @warn "These parameters converged to (or very near) their bounds: $pinned. This usually " *
              "means the fit is under-constrained or the bounds are too tight for this data."
    end
    condition_number
end
