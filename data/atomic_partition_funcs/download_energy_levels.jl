# this file is separate from the script which calculates and saves partition functions because 
# it's sometimes handy to include this file to download energy levels for ad-hoc analysis 

# the is a more sophisticated parser for NIST energy level data in data/metal_bf_cross-sections

using Korg, Downloads, CSV, DataFrames

function download_levels_from_NIST(elems=Korg.atomic_symbols)
    dfs = Dict()
    for spec in Korg.all_atomic_species()
        if spec in Korg.Species.(["H II", "He III"])
            continue
        end

        species = replace(string(spec), ' ' => '+')
        url = "https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&spectrum=$(species)&units=1&format=2&output=0&page_size=15&multiplet_ordered=1&conf_out=on&level_out=on&g_out=on&temp=&submit=Retrieve+Data"
        df = CSV.read(Downloads.download(url), DataFrame;
                      header=["configuration", "g", "level", "c3"], skipto=2)

        df.level = map(df.level) do l
            #drop characters which supply extra info about line
            parse(Float64,
                  replace(l, r"\+." => "", "&dagger;" => "", (collect("\"[]()=au?") .=> "")...))
        end
        select!(df, ["configuration", "g", "level"])
        filter!(df) do row #we don't want any "summary" rows
            !isnothing(match(r"[a-z]", row.configuration))
        end
        dropmissing!(df)

        df.n = map(df.configuration) do conf
            #last instance of an integer after either " or .
            parse(Int, (collect(eachmatch(Regex("[\"\\.](\\d+)"), conf))[end].captures[1]))
        end

        dfs[spec] = df
    end
    dfs
end
