module Fit
using ..Korg, ProgressMeter, ForwardDiff
using SparseArrays: spzeros

"""
TODO
"""
function downsampled_LSF_matrix(synth_wls, obs_wls, R)
    #ideas - require wls to be a rangbe object? Use erf to account for grid edges?
    convM = spzeros((length(obs_wls), length(synth_wls)))
    lb, ub = 1,1 #initialize window bounds
    @showprogress  "Constructing LSF matrix" for i in 1:length(obs_wls)
        λ0 = obs_wls[i]
        σ = λ0 / R / 2
        lb, ub = Korg.move_bounds(synth_wls, lb, ub, λ0, 3σ)
        ϕ = Korg.normal_pdf.(synth_wls[lb:ub] .- λ0, σ) * step(synth_wls)
        @. convM[i, lb:ub] += ϕ
    end
    convM 
end

# scaling free params.
function scale(p)
    nodes = Korg.get_atmosphere_archive()[1]
    lower = first.(nodes)
    upper = last.(nodes)

    p = (p - lower) ./ (upper - lower)
    tan.(π .* (p.-0.5))
end

function unscale(p)
    nodes = Korg.get_atmosphere_archive()[1]
    lower = first.(nodes)
    upper = last.(nodes)

    p = atan.(p)./π .+ 0.5
    p.*(upper - lower) .+ lower
end

function synth(synthesis_wls, obs_wls, scaled_p, linelist, LSF_matrix)
    p = unscale(scaled_p)
    atm = Korg.interpolate_marcs(p...)
    alpha_els = ["Ne", "Mg", "Si", "S", "Ar", "Ca", "Ti"]
    A_X = Korg.format_A_X(p[3], Dict(["C"=>p[5]+p[3], (alpha_els .=> p[4]+p[3])...]))
    F = Korg.synthesize(atm, linelist, A_X, synthesis_wls; vmic=1.0).flux
    #rF = apply_rotation(F, wls, 10) #TODO
    dF = LSF_matrix * F
    Korg.rectify(dF, obs_wls; wl_step=0)
end

"""
TODO rotation, vmic, abundances
"""
function find_best_params(obs_wls, data, err, R, linelist; 
                          p0=(upper+lower)./2, synthesis_wls=first(obs_wls):0.01:last(obs_wls),
                          LSF_matrix=downsampled_LSF_matrix(synthesis_wls, obs_wls, R)
                          )
    
    global_synth = let synthesis_wls=synthesis_wls, obs_wls=obs_wls, linelist=linelist, LSF_matrix=LSF_matrix
        p -> synth(synthesis_wls, obs_wls, p, linelist, LSF_matrix)
    end
    #global_chi2 = let wls=wls, data=data, err=err
    #    p -> begin
    #        flux = synth(wls, unscale(p))
    #        sum(((flux .- data)./err).^2)
    #    end
    #end 
    
    @time J0 = ForwardDiff.jacobian(global_synth, scale(p0))
    sum(J0.^2, dims=2)
end

end # module
