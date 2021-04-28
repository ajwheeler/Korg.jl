push!(LOAD_PATH,"../src/")

using Documenter, SSSynth

makedocs(sitename="SSSynth", 
         modules=[SSSynth],
         pages=[
                "Quickstart" => "index.md"
                "All Functions" => "API.md"
                "References" => "refs.md"
               ],
        authors="Adam Wheeler and Matthew Abruzzo",
        format=Documenter.HTML(assets=["assets/favicon.ico"])
       )

deploydocs(repo = "github.com/ajwheeler/SSSynth.jl.git",
           devbranch="main")
