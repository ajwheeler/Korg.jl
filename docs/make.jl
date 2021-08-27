push!(LOAD_PATH,"../src/")

using Documenter, Korg

makedocs(sitename="Korg", 
         modules=[Korg],
         pages=[
                "Quickstart" => "index.md"
                "Function Reference" => "API.md"
                "Developer Documentation" => "devdocs.md"
                "Install" => "install.md"
               ],
        authors="Adam Wheeler and Matthew Abruzzo",
        format=Documenter.HTML(assets=["assets/favicon.ico"])
       )

deploydocs(repo = "github.com/ajwheeler/Korg.jl.git",
           devbranch="main")
