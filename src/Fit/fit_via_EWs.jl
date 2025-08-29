using DiffResults, Trapz
using Statistics: mean, std
using Optim

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
    uranium ([`Korg.format_A_X`](@ref) exists as a convenience for creating this vector). All
    syntheses are done with these abundances, so if the resulting abundances deviate significantly
    from these, you may wish to iterate.

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
                                typeof(linelist[1].log_gf)) # this is a bit of a hack
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
    ews_to_abundances(atm, linelist, A_X, measured_EWs; kwargs...)
    ews_to_abundances(params, linelist, measured_EWs; kwargs...)

Compute per-line abundances on the linear part of the curve of growth given a model atmosphere and a
list of lines with equivalent widths.

# Arguments:

  - `atm`: the model atmosphere (see [`Korg.read_model_atmosphere`](@ref) and
    [`Korg.interpolate_marcs`](@ref)).
  - `linelist`: A vector of [`Korg.Line`](@ref)s (see [`Korg.read_linelist`](@ref)).  The lines must
    be sorted by wavelength.
  - `A_X`: a vector containing the A(X) abundances (log(n_X/n_H) + 12) for elements from hydrogen to
    uranium (see [`Korg.format_A_X`](@ref)). All syntheses are done with these abundances, so if the
    resulting abundances deviate significantly from these, you may wish to iterate.
  - `measured_EWs`: a vector of equivalent widths (in mÅ)

Alternatively, you can pass a vector of parameters like the one returned by
[`ews_to_stellar_parameters`](@ref) instead of a model atmosphere and abundances vector. These
should be in the order `[Teff, logg, vmic, M_H]`.

# Returns

A vector of abundances (`A(X) = log10(n_X/n_H) + 12` format) for each line in`linelist`, and a
vector of ∂A/∂log(EW) for each line.

# Optional arguments:

  - `wl_step` (default: 0.01) is the resolution in Å at which to synthesize the spectrum around each
    line.

  - `ew_window_size` (default: 2): the farthest (in Å) to consider equivalent width contributions for
    each line.  It's very important that this is large enough to include each line entirely.
  - `blend_warn_threshold` (default: 0.01) is the minimum absorption between two lines allowed before
    triggering a warning that they may be blended.
  - `finite_difference_delta_A` (default: 0.01): the step size in A(X) to use for the finite
    difference calculation of the curve of growth slope.
  - `solar_abundances` (default: `Korg.default_solar_abundances`): the solar abundances to use for
    the calculation of A(X). Only used if `ews_to_abundances` is called with a vector of parameters
    instead of a model atmosphere and abundances vector.

All other keyword arguments are passed to [`Korg.synthesize`](@ref) when synthesizing each line.
"""
function ews_to_abundances(atm, linelist, A_X, measured_EWs; ew_window_size::Real=2.0,
                           wl_step=0.01, callback=Returns(nothing), abundance_tol=1e-5,
                           max_iter=30, blend_warn_threshold=0.01, verbose=false,
                           solar_abundances=Korg.default_solar_abundances,
                           synthesize_kwargs...)
    ews_to_abundances_parameter_validation(linelist, measured_EWs)

    # do a single synthesis to get the chemical equilibrium once
    sol = synthesize(atm, [], A_X, (5000, 5000))

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
    ∂A_∂logEW = similar(ΔA)
    while sum(fitmask) > 0 && iter < max_iter # while there are still lines to fit
        iter += 1

        # calculate the EWs for ΔA[fitmask] and it's derivative wrt log(EW)
        perturbed_linelist = [Korg.Line(l; log_gf=l.log_gf .+ Δ)
                              for (l, Δ) in zip(linelist[fitmask], ΔA[fitmask])]
        EWs[fitmask] .= calculate_EWs(atm, perturbed_linelist, A_X; ew_window_size=ew_window_size,
                                      wl_step=wl_step, blend_warn_threshold=blend_warn_threshold,
                                      use_chemical_equilibrium_from=sol, synthesize_kwargs...)
        ∂A_∂logEW[fitmask] .= @. (ΔA[fitmask] - ΔA_prev[fitmask]) /
                                 log10(EWs[fitmask] / EWs_prev[fitmask])

        # calculate the change in abundances according to the secant method
        δA = @. ∂A_∂logEW[fitmask] * log10(measured_EWs[fitmask] / EWs[fitmask])

        # update the "previous" abundances and EWs
        ΔA_prev[fitmask] .= ΔA[fitmask]
        EWs_prev[fitmask] .= EWs[fitmask]

        # update the abundances of the lines that are still being fit
        ΔA[fitmask] .+= δA

        callback(A0 .+ ΔA)

        # stop fitting lines that have converged
        fitmask[fitmask] .&= abs.(δA) .> abundance_tol

        verbose && println("iter $iter ($(sum(fitmask)) lines unconverged)")
    end

    # NaN out anything unconverged
    ΔA[fitmask] .= NaN

    A0 .+ ΔA
end
function ews_to_abundances(params, linelist, measured_EWs;
                           solar_abundaces=Korg.default_solar_abundances, kwargs...)
    Teff, logg, vmic, M_H = params
    A_X = Korg.format_A_X(M_H; solar_abundances=solar_abundaces)
    atm = Korg.interpolate_marcs(Teff, logg, A_X; perturb_at_grid_values=true,
                                 clamp_abundances=true)
    ews_to_abundances(atm, linelist, A_X, measured_EWs; vmic=vmic, kwargs...)
end

"""
    ews_to_abundances_approx(atm, linelist, A_X, measured_EWs; kwargs...)

A very approximate method for fitting abundances from equivalent widths.  It assumes that all lines
are on the linear part of the curve of growth, and that the lines are not blended. Arguments and
keyword arguments are the same as for [`ews_to_abundances`](@ref). Returns the abundances, but not
∂A/∂logEW, since this is constant by assumption.
"""
function ews_to_abundances_approx(atm, linelist, A_X, measured_EWs; ew_window_size::Real=2.0,
                                  wl_step=0.01, blend_warn_threshold=0.01, synthesize_kwargs...)
    A0 = [A_X[Korg.get_atoms(l.species)[1]] for l in linelist]

    EWs = calculate_EWs(atm, linelist, A_X; ew_window_size=ew_window_size, wl_step=wl_step,
                        blend_warn_threshold=blend_warn_threshold, synthesize_kwargs...)

    @. A0 + (log10(measured_EWs) - log10.(EWs))
end
function ews_to_abundances_approx(params, linelist, measured_EWs;
                                  solar_abundances=Korg.default_solar_abundances, kwargs...)
    Teff, logg, vmic, M_H = params
    A_X = Korg.format_A_X(M_H; solar_abundances=solar_abundances)
    atm = Korg.interpolate_marcs(Teff, logg, A_X; perturb_at_grid_values=true,
                                 clamp_abundances=true)
    ews_to_abundances_approx(atm, linelist, A_X, measured_EWs; vmic=vmic, kwargs...)
end

# basic checks for the parameters of ews_to_abundances and ews_to_abundances_approx
function ews_to_abundances_parameter_validation(linelist, measured_EWs)
    if length(linelist) != length(measured_EWs)
        throw(ArgumentError("Length of linelist does not match length of measured EWs ($(length(linelist)) != $(length(measured_EWs)))"))
    end

    if any(l -> Korg.ismolecule(l.species), linelist)
        throw(ArgumentError("Linelist contains molecular species. Only atomic lines can be used when fitting EWs."))
    end

    # Check that the user is supplying EWs in mA
    if 1 > maximum(measured_EWs)
        @warn "Maximum EW given is less than 1 mA. Check that you're giving EWs in mÅ (*not* Å)."
    end
end

"""
    ews_to_stellar_parameters_direct(linelist, measured_EWs,
                                    measured_EW_err=ones(length(measured_EWs));
                                    p0=[5.0, 3.5, 1.0, 0.0], precision=1e-5, time_limit=500)

!!! warning

    This function is experimental and may change or be removed in the future.

Find stellar parameters from equivalent widths the "direct" way, e.g. by forward modelling the
equivalent widths and minimizing the chi-squared.

Returns a vector of parameters `[Teff, logg, vmic, [m/H]]`, and an estimate of the uncertainties
from the approximate inverse Hessian matrix from BFGS.
"""
function ews_to_stellar_parameters_direct(linelist, measured_EWs,
                                          measured_EW_err=ones(length(measured_EWs));
                                          verbose=false, p0=[5.0, 3.5, 1.0, 0.0], precision=1e-5,
                                          time_limit=500)
    function cost(params)
        Teff, logg, vmic, M_H = params
        Teff = Teff * 1e3
        A_X = format_A_X(M_H)
        atm = interpolate_marcs(Teff, logg, A_X)
        EWs = calculate_EWs(atm, linelist, A_X; vmic=vmic)
        sum(@. ((EWs - measured_EWs) / measured_EW_err)^2)
    end

    @time res = optimize(cost, p0, BFGS(; linesearch=LineSearches.BackTracking(; maxstep=1.0)),
                         Optim.Options(; x_abstol=precision, time_limit=time_limit,
                                       extended_trace=true, store_trace=true, show_trace=verbose);
                         autodiff=:forward)

    if !Optim.converged(res)
        @warn "Stellar parameter fit to EWs did not converge"
    end
    params = res.minimizer
    params[1] *= 1e3

    scales = [1e3, 1.0, 1.0, 1.0]
    uncertainties = res.trace[end].metadata["~inv(H)"]
    uncertainties .*= scales .* scales'

    params, uncertainties
end

function _stellar_params_default_callback(params, residuals, abundances)
    println("params: ", params)
    println("residuals: ", residuals)
    println()
end

"""
    ews_to_stellar_parameters(linelist, measured_EWs; kwargs...)

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

# Returns:

A pair, `(params, uncertainties)`, containing:

  - the best-fit parameters: `[Teff, logg, vmic, [m/H]]` as a vector (with vmic in km/s)
  - the uncertainties in the parameters, propagated from the uncertainties in the
    equivalent widths, which estimated from the line-to-line scatter and assumed to be uncorrelated.
    This is not a particularly rigorous error estimate, but it is a good sanity check.

!!! info

    For now, this function is limited to the parameters supported by the SDSS MARCS grid, i.e.
    down to [m/H] = -2.5 at the lowest metallicity grid point.

# Keyword arguments:

  - `abundance_adjustments` (default: `zeros(length(measured_EWs))`) is a vector of abundances
    adjustments to be applied to the lines. This is useful when doing a differential analysis.

  - `Teff0` (default: 5000.0) is the starting guess for Teff
  - `logg0` (default: 3.5) is the starting guess for logg
  - `vmic0` (default: 1.0) is the starting guess for vmic. Note that this must be nonzero in order to
    avoid null derivatives. Very small values are fine.
  - `M_H0` (default: 0.0) is the starting guess for [m/H]
  - `tolerances` (default: `[1e-3, 1e-3, 1e-4, 1e-3]`) is the tolerance for the residuals each equation
    listed above. The solver stops when all residuals are less than the corresponding tolerance.
  - `max_step_sizes` (default: `[1000.0, 1.0, 0.3, 0.5]`) is the maximum step size to take in each
    parameter direction.  This is used to prevent the solver from taking too large of a step and
    missing the solution.  Be particularly cautious with the vmic (third) parameter, as the
    unadjusted step size is often too large.
  - `parameter_ranges` (default: `[(2800.0, 8000.0), (-0.5, 5.5), (1e-3, 10.0), (-2.5, 1.0)]`) is the
    allowed range for each parameter. This is used to prevent the solver from wandering into
    unphysical parameter space, or outside the range of the SDSS MARCS grid supported by
    [`Korg.interpolate_marcs`](@ref). The default ranges for ``T_\\mathrm{eff}``, ``\\log g``, and [m/H]
    are the widest supported by the MARCS grid. The minimum value for `vmic` is set to 1e-3 km/s to
    avoid null derivatives in the optimization.
  - `fix_params` (default: `[false, false, false, false]`) is a vector of booleans indicating which
    parameters should be held fixed during the optimization. The order is [Teff, logg, vmic, [m/H]].
  - `solar_abundances` (default: [`Korg.default_solar_abundances`](@ref)) the solar abundances to assume
  - `verbose` (default: `false`) is a boolean indicating whether to print verbose output during the
    optimization. Works, by default, by adding a callback that prints the current parameters,
    residuals, and abundances. The callback can be explicitly set with the `callback` keyword.
  - `callback`: is a function which is called at each step of the optimization.
    It is passed three arguments:

      + the current values of the parameters [Teff, logg, vmic, [m/H]]
      + the residuals of each equation being solved
      + the abundances of each line computed with the current parameters
        You can pass a callback function, to e.g. make a plot of the residuals at each step.
  - `max_iterations` (default: 30) is the maximum number of iterations to allow before stopping the
    optimization.
"""
function ews_to_stellar_parameters(linelist, measured_EWs;
                                   abundance_adjustments=zeros(length(measured_EWs)),
                                   Teff0=5000.0, logg0=3.5, vmic0=1.0, M_H0=0.0,
                                   tolerances=[1e-3, 1e-3, 1e-4, 1e-3],
                                   max_step_sizes=[1000.0, 1.0, 0.3, 0.5],
                                   parameter_ranges=[(2800.0, 8000.0),
                                       (-0.5, 5.5),
                                       (1e-3, 10.0),
                                       (-2.5, 1.0)],
                                   fix_params=[false, false, false, false],
                                   solar_abundances=Korg.default_solar_abundances,
                                   verbose=false,
                                   callback=if verbose
                                       _stellar_params_default_callback
                                   else
                                       Returns(nothing)
                                   end,
                                   max_iterations=30, passed_kwargs...)
    if :vmic in keys(passed_kwargs)
        throw(ArgumentError("vmic must not be specified, because it is a parameter fit by " *
                            "ews_to_stellar_parameters. Did you mean to specify vmic0, the starting " *
                            "value? See the documentation for ews_to_stellar_parameters if you would " *
                            "like to fix microturbulence to a given value."))
    end

    params0 = [Teff0, logg0, vmic0, M_H0]
    params = validate_ews_to_stellar_params_inputs(linelist, measured_EWs, params0,
                                                   parameter_ranges)

    # First phase: approximate calculations
    get_residuals = (p) -> _stellar_param_equation_residuals(false, p, linelist, measured_EWs,
                                                             abundance_adjustments,
                                                             solar_abundances, fix_params,
                                                             callback, passed_kwargs)
    J_result = DiffResults.JacobianResult(params)
    if !_ews_to_stellar_parameters_phase!(J_result, get_residuals, params, fix_params,
                                          tolerances * 10, max_step_sizes, parameter_ranges,
                                          max_iterations)
        return params, fill(NaN, 4)
    end

    verbose && println("Approximate solve converged. Starting exact solve.")

    # Second phase: exact calculations
    get_residuals = (p) -> _stellar_param_equation_residuals(true, p, linelist, measured_EWs,
                                                             abundance_adjustments,
                                                             solar_abundances, fix_params,
                                                             callback, passed_kwargs)
    if !_ews_to_stellar_parameters_phase!(J_result, get_residuals, params, fix_params,
                                          tolerances, max_step_sizes, parameter_ranges,
                                          max_iterations)
        return params, fill(NaN, 4)
    end

    # compute uncertainties
    residual_uncertainties = _stellar_param_residual_uncertainties(params, linelist, measured_EWs,
                                                                   abundance_adjustments,
                                                                   solar_abundances, passed_kwargs)

    J = DiffResults.jacobian(J_result)[.!fix_params, .!fix_params]
    param_uncertainties = zeros(4)
    param_uncertainties[.!fix_params] .= abs.(J \ residual_uncertainties[.!fix_params])

    params, param_uncertainties
end

# Check lots of things and throw informative errors.
# Clamp the parameters to the specified ranges and warn if necessary.
function validate_ews_to_stellar_params_inputs(linelist, measured_EWs, params0,
                                               parameter_ranges)
    if length(linelist) != length(measured_EWs)
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

    params
end

# Perform one phase of the stellar parameter iteration, either with approximate or exact calculations.
# Returns true if the iteration converged, false otherwise.
function _ews_to_stellar_parameters_phase!(J_result, get_residuals, params, fix_params,
                                           tolerances, max_step_sizes, parameter_ranges,
                                           max_iterations, damping_factor=1.0)
    iterations = 0
    while true
        J_result = ForwardDiff.jacobian!(J_result, get_residuals, params)
        J = DiffResults.jacobian(J_result)
        residuals = DiffResults.value(J_result)
        if all((abs.(residuals).<tolerances)[.!fix_params]) # stopping condition
            return true
        end

        # calculate step, and update params
        step = zeros(length(params))
        step[.!fix_params] = -J[.!fix_params, .!fix_params] \ residuals[.!fix_params]
        step .*= damping_factor
        params .+= clamp.(step, -max_step_sizes, max_step_sizes)
        params .= clamp.(params, first.(parameter_ranges), last.(parameter_ranges))

        iterations += 1
        if iterations > max_iterations
            @warn "Failed to converge after $max_iterations iterations.  Returning the current guess."
            return false
        end
    end
end

# the residuals for the four excitation-ionization equations, using either exact or approximate
# calculations.
function _stellar_param_equation_residuals(exact_calculation, params, linelist, EW,
                                           abundance_adjustments, solar_abundances, fix_params,
                                           callback, passed_kwargs)
    A, neutrals, REWs, Z = _stellar_param_equations_precalculation(exact_calculation, params,
                                                                   linelist, EW,
                                                                   abundance_adjustments,
                                                                   solar_abundances,
                                                                   passed_kwargs)
    finitemask = isfinite.(A)

    if mean(finitemask) < 0.7
        throw(ErrorException("Less than 70% of the lines converged.  This may be due to a " *
                             "problem with the linelist, the measured EWs, or the initial parameter " *
                             "guess."))
    end
    if mean(finitemask[.!neutrals]) < 0.5
        throw(ErrorException("Less than 50% of the ion lines converged.  This may be due to a " *
                             "problem with the linelist, the measured EWs, or the initial parameter " *
                             "guess."))
    end

    neutrals = neutrals .& finitemask

    teff_residual = get_slope([line.E_lower for line in linelist[neutrals]], A[neutrals])
    logg_residual = mean(A[neutrals]) - mean(A[.!neutrals.&finitemask])

    vmic_residual = get_slope(REWs[neutrals], A[neutrals])
    feh_residual = mean(A[finitemask]) - (params[4] + solar_abundances[Z])
    residuals = [teff_residual, logg_residual, vmic_residual, feh_residual]
    residuals .*= .!fix_params # zero out residuals for fixed parameters

    callback(ForwardDiff.value.(params), ForwardDiff.value.(residuals), ForwardDiff.value.(A))
    residuals
end

# called by _stellar_param_equation_residuals
# returns (statistical_uncertainty, systematic_uncertainty)
function _stellar_param_residual_uncertainties(params, linelist, EW, abundance_adjustments,
                                               solar_abundances, passed_kwargs)
    A, neutrals, REWs, _ = _stellar_param_equations_precalculation(true, params, linelist, EW,
                                                                   abundance_adjustments,
                                                                   solar_abundances,
                                                                   passed_kwargs)
    # assume the total (including systematic) err in the abundances of each line can be obtained
    # from the line-to-line scatter
    estimated_err = std(A[isfinite.(A)])

    sigma_mean = estimated_err ./ sqrt(length(estimated_err))
    teff_residual_sigma = estimated_err *
                          get_slope_uncertainty([line.E_lower for line in linelist[neutrals]])
    vmic_residual_sigma = estimated_err * get_slope_uncertainty(REWs[neutrals])

    [teff_residual_sigma, sigma_mean, vmic_residual_sigma, sigma_mean]
end

function _stellar_param_equations_precalculation(exact_calculation, params, linelist, EW,
                                                 abundance_adjustments, solar_abundances,
                                                 passed_kwargs)
    if exact_calculation
        A = ews_to_abundances(params, linelist, EW; blend_warn_threshold=Inf,
                              solar_abundances=solar_abundances, passed_kwargs...)
    else
        A = ews_to_abundances_approx(params, linelist, EW; blend_warn_threshold=Inf,
                                     solar_abundances=solar_abundances, passed_kwargs...)
    end

    neutrals = [l.species.charge == 0 for l in linelist]
    REWs = log10.(EW ./ [line.wl for line in linelist])
    # this is guaranteed not to be a mol (checked by ews_to_stellar_parameters).
    Z = Korg.get_atoms(linelist[1].species)[1]

    A + abundance_adjustments, neutrals, REWs, Z
end

# called by _stellar_param_equation_residuals
function get_slope(xs, ys)
    Δx = xs .- mean(xs)
    Δy = ys .- mean(ys)
    sum(Δx .* Δy) ./ sum(Δx .^ 2)
end
function get_slope_uncertainty(xs)
    n = length(xs)
    sqrt(1 / (sum(xs .^ 2) - sum(xs)^2 / n))
end
