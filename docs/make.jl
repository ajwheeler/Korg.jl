using Documenter, Korg

# Use the README.md as index.md in the docs. These should really have the same contents.
cp("../README.md", "./src/index.md"; force=true)

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

         # ideally we would do link checking, but adam doesn't like it turning his badge red
         linkcheck=false, # check if any external links are broken
         # there are useful if linkchecking gets turned back on
         linkcheck_timeout=30, # default (10s) is too short for slow academic sites
         linkcheck_ignore=["https://marcs.astro.uu.se/"],

         authors="Adam Wheeler, Matthew Abruzzo, Andrew Casey, and collaborators",
         format=Documenter.HTML(; assets=["assets/favicon.ico"], size_threshold=500000),
         checkdocs=:all,)

deploydocs(; repo="github.com/ajwheeler/Korg.jl.git", devbranch="main")
