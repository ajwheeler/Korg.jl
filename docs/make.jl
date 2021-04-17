push!(LOAD_PATH,"../src/")

using Documenter, SSSynth

makedocs(sitename="SSSynth Documentation", 
         modules=[SSSynth],
         pages=[
                "Quickstart" => "index.md",
                "All Functions" => "API.md"
               ])

deploydocs(repo = "github.com/ajwheeler/SSSynth.jl.git",
           devbranch="main")
