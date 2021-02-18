include("constants.jl")

blackbody(λ, T) = 2*h*c^2/λ^5 * 1/(exp(h*c/λ/k/T) - 1)


"""trapezoid rule"""
function trap(xs, f)
    sum((f[1:end-1] + 0.5diff(f)) .* diff(xs))
end

"""
solve transport equation
"""
function intensity(λs, falc, α=0, μ=1)

    αHmin(λ) = @. extinction_Hmin(λ, falc.temp, falc.n_electron * k * falc.temp) * (falc.n_Htotal - falc.n_proton)

    map(λs) do λ
        α = extinction_thompson.(falc.n_electron) + αHmin(λ)
        τ = zeros(size(falc.height))
        for i in 2:length(τ)
            τ[i] = τ[i-1] + 0.5(α[i] + α[i-1])*(falc.height[i-1] - falc.height[i])
        end
        trap(τ./μ, blackbody.(λ, falc.temp).*exp.(-τ/μ))
    end
end

"""
compute astrophysical flux
"""
function astrophysical_flux(z, α; ϵ=1e-6)
    #TODO combine with depth integration?
    ans, = quadgk(μ->intensity(λ, μ), ϵ, 1)
    ans
end

"extinction from thompson scattering"
extinction_thompson(nₑ)  = 6.648e-25 * nₑ

"""
out: H-minus bf+ff extinction [cm^2 per neutral hydrogen atom]
     assuming LTE ionization H/H-min
     NB: includes stimulated emission correction already!
"""
function extinction_Hmin(λ,T,pₑ)
    θ = 5040. / T
    wav = λ*1e8 #work in aangstroms

    # evaluate H-min bound-free per H-min ion = Gray (8.11)
    # his alpha = my sigma in NGSB/AFYC (per particle without stimulated)
    if wav > 16444
        σbf = 0.
    else
        σbf = 1.99654 - 1.18267e-5*wav + 2.64243e-6*wav^2 - 4.40524e-10*wav^3 +
                3.23992e-14*wav^4 -1.39568e-18*wav^5 +2.78701e-23*wav^6
        σbf *= 1e-18  # cm^2 per H-min ion
    end

    # convert into bound-free per neutral H atom assuming Saha = Gray p135
    # units: cm2 per neutral H atom in whatever level (whole stage)
    graysaha = 4.158e-10 * pₑ * θ^2.5 * 10^(0.754θ) # Gray (8.12)
    κbf=σbf*graysaha                               # per neutral H atom
    κbf *= (1 - exp(-h*c/(wav*1e-8*k*T)))
    
    # evaluate H-min free-free including stimulated emission = Gray p136
    lwav = log10(wav)
    f0 =  -2.2763 -1.6850*lwav +0.76661*lwav^2 -0.0533464*lwav^3
    f1 =  15.2827 -9.2846*lwav +1.99381*lwav^2 -0.142631*lwav^3
    f2 = -197.789 +190.266*lwav -67.9775*lwav^2 +10.6913*lwav^3 -0.625151*lwav^4                                                    
    ltheta = log10(θ)
    κff = 1e-26 * pₑ * 10^(f0+f1*ltheta+f2*ltheta^2)   # Gray (8.13)

    κbf + κff #TODO - not clear why these are κs and not αs
end

