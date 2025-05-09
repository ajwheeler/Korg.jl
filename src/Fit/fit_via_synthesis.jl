using LineSearches, Optim
using Interpolations: linear_interpolation, Line

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
    # wavelengths, windows and LSF
    synthesis_wls, obs_wl_mask, LSF_matrix = _setup_wavelengths_and_LSF(obs_wls, synthesis_wls,
                                                                        LSF_matrix, R, windows,
                                                                        wl_buffer)

    _validate_observed_spectrum(obs_wls, obs_flux, obs_err, obs_wl_mask)

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
