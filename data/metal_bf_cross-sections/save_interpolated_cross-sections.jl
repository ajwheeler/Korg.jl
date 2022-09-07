using Korg, HDF5, Base
using Statistics: mean

include("NORAD_TOPBase_parsing_code.jl")

DATA_DIR = ARGS[1]
OUT_FILE = ARGS[2]

function save_bf_cross_section(f, spec, logTs, νs, ξ=2e5)
    @assert issorted(νs, rev=true) # should be decreasing freq
    
    Ts = 10 .^ logTs
    wls = Korg.c_cgs ./ νs
    σs = cross_section_bf(spec, wls, Ts, DATA_DIR)
    
    Rs = Korg.c_cgs ./ sqrt.(2 * Korg.kboltz_cgs * Ts ./ Korg.get_mass(Korg.species"He I") .+ ξ^2) 
    λs = Korg.c_cgs ./ νs # the frequencecy grid is linear-enough in wavelength over the ~ 1 Å scale LSF 

    # convert to megabarns and cast to single precision
    σs = hcat(Korg.constant_R_LSF.(eachcol(σs), Ref(λs), Rs)...)
    mask = σs .< floatmin(Float32)
    σs[mask] .= 0
    print(spec)
    println(" ", minimum(σs), " -- ", maximum(σs), " ($(mean(mask)) denormal)")
    log_σs = Float32.(log.((σs)))
    
    # save in order of increasing freq / decreasing wl
    f["cross-sections/"*string(spec), shuffle=(), deflate=6] = reverse(log_σs, dims=1)
end

logTs = 2:0.1:5
wl_lo = 490 * 1e-8
wl_hi = 30_050 * 1e-8

# descending freq, ascending wl
νs = Korg.c_cgs / wl_lo : -1e11 : Korg.c_cgs / wl_hi

filename = OUT_FILE
h5open(filename, "w") do f
    f["logT_min"]  = logTs[1]
    f["logT_step"] = step(logTs)
    f["logT_max"] = logTs[end]
    f["nu_min"] = νs[end]
    f["nu_step"] = -step(νs)
    f["nu_max"] = νs[1]

    atomic_numbers = [1:14... ; 16:2:20... ; 26]
    @time for Z in atomic_numbers
        for ionization_number in [1, 2]
            if Z == 1 && ionization_number == 2
                continue # no such thing as H II 
            end
            spec = Korg.Species(Korg.atomic_symbols[Z] * " $ionization_number")
            save_bf_cross_section(f, spec, logTs, νs)
        end
    end
end
;
