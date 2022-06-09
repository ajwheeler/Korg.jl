# this file is separate from the script which calculates and saves partition functions because 
# it's sometimes handy to include this file to download energy levels for ad-hoc analysis 

using Korg, Downloads, CSV, DataFrames

function download_levels_from_NIST()
    dfs = Dict()
    for elem in Korg.atomic_symbols, ionization in ["I", "II", "III"]
        if (elem == "H" && ionization != "I") || (elem == "He" && ionization == "III")
            continue
        end
        
        species = elem*"+"*ionization
        
        url = "https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&spectrum=$(species)&units=1&format=2&output=0&page_size=15&multiplet_ordered=1&conf_out=on&level_out=on&g_out=on&temp=&submit=Retrieve+Data"
        df = CSV.read(Downloads.download(url), DataFrame, header=["configuration", "g", "level", "c3"], skipto=2)
        
        df.level = map(df.level) do l
            #drop characters which supply extra info about line
            parse(Float64, replace(l, r"\+."=>"", "&dagger;"=>"", (collect("\"[]()=au?") .=> "")...))
        end
        select!(df, ["configuration", "g", "level"])
        dropmissing!(df)
        df = df[df.configuration .!= "=\"\"", :]
        
        df.n = map(df.configuration) do conf
            #last instance of an integer after either " or .
            parse(Int, (collect(eachmatch(Regex("[\"\\.](\\d+)"), conf))[end].captures[1]))
        end
        
        dfs[Korg.Species("$elem $ionization")] = df
    end
    dfs
end