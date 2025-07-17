# this script produces solar abundances from Bergemann et al. (2025) in a form that
# can be copy-pasted src/atomic_data.jl, which I hope is grug-brained is a good way

using CSV, DataFrames
df = CSV.read("Table1_Photosphere_Bergemann_Lodders_Palme_2025_Abundances.csv", DataFrame; header=2,
              skipto=5, footerskip=1, silencewarnings=true)
rename!(df, " Sun Convection Zone" => :sun, " CI-Chondrites" => :chondrite)

for col in [:sun, :chondrite]
    df[!, col] = map(df[!, col]) do element
        if element == "..."
            missing
        else
            parse(Float64, element)
        end
    end
end

abunds = Dict(df.Z .=> tuple.(df.sun, df.chondrite))

map(1:92) do Z
    if Z in keys(abunds)
        A_conv, A_chond = abunds[Z]
        if !ismissing(A_conv)
            A_conv
        elseif !ismissing(A_chond)
            A_chond
        else
            -5.0
        end
    else
        -5.0
    end
end |> println
