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
         checkdocs=:all, # make sure all docstrings are in the documentation
         linkcheck=true, # check if any external links are broken 
         linkcheck_timeout=30, # default (10s) is too short for slow academic sites
         linkcheck_ignore=[], # maybe don't check ADS links because they like to time out
         authors="Adam Wheeler and Matthew Abruzzo",
         format=Documenter.HTML(; assets=["assets/favicon.ico"]))

deploydocs(; repo="github.com/ajwheeler/Korg.jl.git",
           devbranch="main")
