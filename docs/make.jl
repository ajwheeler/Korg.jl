#push!(LOAD_PATH, "../src/")

using Documenter, Korg

makedocs(;
         modules=[Korg],
         repo=Documenter.Remotes.GitHub("ajwheeler", "Korg.jl"),
         sitename="Korg",
         pages=["Quickstart" => "index.md"
                "Install" => "install.md"
                "Tutorials" => "tutorials.md"
                "High-level API" => "API.md"
                "FAQ" => "FAQ.md"
                "Changelog" => "changelog.md"
                "Developer Documentation" => "devdocs.md"],
         linkcheck=true, # check if any external links are broken 
         linkcheck_timeout=30, # default (10s) is too short for slow academic sites
         # dodgy academic links sometimes timeout
         linkcheck_ignore=["https://marcs.astro.uu.se/"],
         authors="Adam Wheeler, Matthew Abruzzo, Andrew Casey, and collaborators",
         format=Documenter.HTML(; assets=["assets/favicon.ico"], size_threshold=500000),)
checkdocs = :all, # make sure all docstrings are in the documentation
            deploydocs(; repo="github.com/ajwheeler/Korg.jl.git",
                       devbranch="main")
