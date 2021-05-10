push!(LOAD_PATH,"../src/")

using Documenter, Korg

makedocs(sitename="Korg", 
         modules=[Korg],
         pages=[
                "Quickstart" => "index.md"
                "Install" => "install.md"
                "All Functions" => "API.md"
                "References" => "refs.md"
               ],
        authors="Adam Wheeler and Matthew Abruzzo",
        format=Documenter.HTML(assets=["assets/favicon.ico"])
       )

deploydocs(repo = "github.com/ajwheeler/Korg.jl.git",
           devbranch="main")
