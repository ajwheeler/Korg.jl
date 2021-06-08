using Random
using Optim

function best_fit_params(atm, linelist, wls, f_obs; ivar=ones(length(f_obs)),
                         free_parameters=["metallicity", "vmic"], fixed_values::Dict=Dict(),
                         initial_values::Dict=Dict(), wl_deflation_factor=0.01)

    function chi2(params, wl_mask)
        f_synth = synthesize(atm, linelist, wls[wl_mask]; metallicity=params[1],
                             abundances=Dict(["Ca"=>params[2]]), verbose=false).flux
        f_synth = Korg.constant_R_LSF(f_synth, wls[wl_mask], 11500.0)
        sum((f_obs[wl_mask] .- f_synth).^2 .* ivar[wl_mask])
    end
    
    mask = rand(length(wls)) .< wl_deflation_factor
    
    result = optimize(p->chi2(p, mask), [0.0, 0.0], NewtonTrustRegion(),
        Optim.Options(show_trace=true, extended_trace=true, x_tol=1e-1);
        autodiff=:forward)
    
    optimize(p->chi2(p, ones(Bool, length(wls))), result.minimizer, NewtonTrustRegion(),
        Optim.Options(show_trace=true, extended_trace=true, x_tol=1e-3);
        autodiff=:forward)
end
