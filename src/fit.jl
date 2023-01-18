module Fit
using ..Korg, ProgressMeter, ForwardDiff
using SparseArrays: spzeros
using Statistics: mean
using LineSearches, Optim
#TODO specify versions in Project.toml

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
    F = try
        Korg.synthesize(atm, linelist, A_X, synthesis_wls; vmic=1.0).flux
    catch e
        ones(length(synthesis_wls))
    end
    #rF = apply_rotation(F, wls, 10) #TODO
    dF = LSF_matrix * F
    Korg.rectify(dF, obs_wls; wl_step=0)
end

"""
TODO rotation, vmic, abundances
"""
function find_best_params(obs_wls, data, err, R, linelist; 
                          p0=(upper+lower)./2, 
                          synthesis_wls=(first(obs_wls)-10) : 0.01 : (last(obs_wls)+10),
                          LSF_matrix=downsampled_LSF_matrix(synthesis_wls, obs_wls, R),
                          verbose=true
                          )
    
    global_synth = let synthesis_wls=synthesis_wls, obs_wls=obs_wls, linelist=linelist, LSF_matrix=LSF_matrix
        p -> synth(synthesis_wls, obs_wls, p, linelist, LSF_matrix)
    end
    chi2 = let obs_wls=obs_wls, data=data, err=err, synthesis_wls=synthesis_wls, LSF_matrix=LSF_matrix, linelist=linelist
        scaled_p -> begin
            flux = synth(synthesis_wls, obs_wls, scaled_p, linelist, LSF_matrix)
            sum(((flux .- data)./err).^2)
        end
    end 

    #@time result = optimize(chi2, scale(p0), BFGS(linesearch=LineSearches.BackTracking()),  
    #                        Optim.Options(x_tol=1e-5, time_limit=3600, store_trace=true, extended_trace=true),
    #                        autodiff=:forward)

    @time J0 = ForwardDiff.jacobian(global_synth, scale(p0))
    fit_each_param_locally(J0, obs_wls, data, err, linelist, p0, synthesis_wls, LSF_matrix, verbose)

    #grad2 = sum(J0.^2, dims=2)
    #wl_chunck_edges = first(obs_wls) : 20 : last(obs_wls)
    #lb, ub = 1, 1
    #mean_grad2 = []
    #for i in eachindex(wl_chunck_edges[1:end-1])
    #    wl_lb = wl_chunck_edges[i]
    #    wl_ub = wl_chunck_edges[i+1]
    #    lb, ub = Korg.move_bounds(obs_wls, lb, ub, (wl_lb+wl_ub)/2, step(wl_chunck_edges))
    #    push!(mean_grad2, (lb, ub, mean(grad2[lb:ub]))) 
    #end
    #lb, ub, _ = mean_grad2[argmax(last.(mean_grad2))]
    #verbose && @info "the subspectrum with the strongest sensitivity to the parameters is $(obs_wls[lb]) – $(obs_wls[ub]) Å"
    #find_best_params_locally(lb, ub, obs_wls, data, err, linelist, p0, synthesis_wls, LSF_matrix, 
    #                         verbose)

    #grad2, lb, ub
end

function fit_each_param_locally(J, obs_wls, data, err, linelist, p0, synthesis_wls, LSF_matrix, verbose; wl_buffer=10)
    wl_chunck_edges = first(obs_wls) : 20 : last(obs_wls)
    lb, ub = 1, 1
    map(eachindex(p0)) do i
        mean_J = []
        for wl_ind in eachindex(wl_chunck_edges[1:end-1])
            wl_lb = wl_chunck_edges[wl_ind]
            wl_ub = wl_chunck_edges[wl_ind+1]
            lb, ub = Korg.move_bounds(obs_wls, lb, ub, (wl_lb+wl_ub)/2, step(wl_chunck_edges))
            push!(mean_J, (lb, ub, mean(J[lb:ub, i])))
        end
        lb, ub, _ = mean_J[argmax(last.(mean_J))]
        verbose && @info "the subspectrum most sensitive to parameter $(i) is $(obs_wls[lb]) – $(obs_wls[ub])"

        synth_wl_lb = findfirst(synthesis_wls .> obs_wls[lb]-wl_buffer)
        synth_wl_ub = findfirst(synthesis_wls .> obs_wls[ub]+wl_buffer) - 1
        small_synthesis_wls = synthesis_wls[synth_wl_lb:synth_wl_ub]
        small_LSF_matrix = LSF_matrix[lb:ub, synth_wl_lb:synth_wl_ub]

        chi2 = let obs_wls=obs_wls[lb:ub], data=data[lb:ub], err=err[lb:ub], p=scale(p0), i=i,
                   synthesis_wls=small_synthesis_wls, LSF_matrix=small_LSF_matrix, linelist=linelist
            scaled_p_i -> begin
                p = [p[1:i-1] ; scaled_p_i ; p[i+1:end]]
                flux = synth(synthesis_wls, obs_wls, p, linelist, LSF_matrix)
                sum(((flux .- data)./err).^2)
            end
        end 

        @time result = optimize(chi2, scale(p0)[i:i], BFGS(linesearch=LineSearches.BackTracking()),  
                                Optim.Options(time_limit=1000, store_trace=true, extended_trace=true))
                                #autodiff=:forward)
    end
end

function find_best_params_locally(lb, ub, obs_wls, data, err, linelist, p0, synthesis_wls, 
                                  LSF_matrix, verbose; wl_buffer=10)
    # synthesize on a range which overshoots the subspectrum by wl_buffer Å on each side
    synth_wl_lb = findfirst(synthesis_wls .> obs_wls[lb]-wl_buffer)
    synth_wl_ub = findfirst(synthesis_wls .> obs_wls[ub]+wl_buffer) - 1
    small_synthesis_wls = synthesis_wls[synth_wl_lb:synth_wl_ub]
    small_LSF_matrix = LSF_matrix[lb:ub, synth_wl_lb:synth_wl_ub]

    chi2 = let obs_wls=obs_wls[lb:ub], data=data[lb:ub], err=err[lb:ub], 
              synthesis_wls=small_synthesis_wls, LSF_matrix=small_LSF_matrix, linelist=linelist
        scaled_p -> begin
            flux = synth(synthesis_wls, obs_wls, scaled_p, linelist, LSF_matrix)
            sum(((flux .- data)./err).^2)
        end
    end 
    verbose && @info "finding best-fit params with subspectrum..."
    # TODO show time if verbose
    @time result = optimize(chi2, scale(p0), BFGS(linesearch=LineSearches.BackTracking()),  
                        Optim.Options(time_limit=1000, store_trace=true, extended_trace=true),
                        autodiff=:forward)
end

end # module
