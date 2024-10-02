push!(LOAD_PATH, "../src/")

using Documenter, Korg

makedocs(; sitename="Korg",
         modules=[Korg],
         pages=["Quickstart" => "index.md"
                "Install" => "install.md"
                "Tutorials" => "tutorials.md"
                "High-level API" => "API.md"
                "FAQ" => "FAQ.md"
                "Changelog" => "changelog.md"
                "Developer Documentation" => "devdocs.md"],
         checkdocs=:all, # make sure all docstrings are in the documentation
         #linkcheck = true, # check if any external links are broken (too many ADS links time out.)
         #linkcheck_timeout = 30, # default (10s) is too short for slow academic sites
         linkcheck_ignore=[],
         strict=[:doctest,
             :linkcheck,
             :parse_error,
             :example_block,
             # This will error on duplicate functions docs, which we currently have by design
             #:autodocs_block, 
             :cross_references,
             :docs_block,
             :eval_block,
             :example_block,
             :footnote,
             :meta_block,
             :missing_docs,
             :setup_block],
         authors="Adam Wheeler and Matthew Abruzzo",
         format=Documenter.HTML(; assets=["assets/favicon.ico"]))

deploydocs(; repo="github.com/ajwheeler/Korg.jl.git",
           devbranch="main")
