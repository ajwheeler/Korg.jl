using CSV, DataFrames

nuclides = CSV.read("BNL_nuclear_wallet_card.dat", DataFrame)

cols = ["A"
        "Element"
        "Z"
        "N"
        "Energy"
        "JPi"
        "Mass Exc"
        "Mass Exc Unc"
        "T1/2 (txt)"
        "T1/2 (seconds)"
        "Abund"
        "Abund_Unc"
        "Dec Mode"
        "Branching (%)"]
rename!(nuclides, cols)

nuclides.JPi = strip.(nuclides.JPi)

#take only nuclides with established spin
filter!(nuclides) do row
    !in(strip(row.JPi), ["", "-"]) &&
        !any(occursin.(['[', 'G', 'g', ',', '(', ')', 'H', 'L'], Ref(row.JPi)))
end

# calculate the numclear spin degeneracy
nuclides.g_ns = map(nuclides.JPi) do JPi
    if in(JPi[end], ['+', '-'])
        JPi = JPi[1:end-1]
    end
    nums = parse.(Int, split(strip(JPi), '/'))
    @assert length(nums) <= 2
    if length(nums) == 1
        (2 * nums[1] + 1)
    else
        @assert nums[2] == 2
        nums[1] + 1
    end
end

# put into dict of dicts
spins = map(pairs(groupby(nuclides, :Z))) do (key, group)
    key.Z => Dict(group.A .=> group.g_ns)
end |> Dict

# pretty print as nicely-formatted Dict 
println("Dict(")
for Z in sort(collect(keys(spins)))
    As = sort(collect(keys(spins[Z])))

    print("    $(Z) => Dict(")
    for A in As[1:end-1]
        print(A, "=>", spins[Z][A], ", ")
    end
    print(As[end], "=>", spins[Z][As[end]])
    if Z == 92
        println(")")
    else
        println("),")
    end
end
println(")")
