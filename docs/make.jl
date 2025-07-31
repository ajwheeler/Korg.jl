using Documenter, Literate, Suppressor, Korg

# Check if we're in docs/ directory or root directory and set this path appropriately
# Locally, I often build from docs, but the CI builds from root
docs_base = basename(pwd()) == "docs" ? "." : "./docs"

# Use the README.md as index.md in the docs.
readme_path = joinpath(docs_base, "..", "README.md")
target_path = joinpath(docs_base, "src", "index.md")
cp(readme_path, target_path; force=true)

# Generate the mardown pages with Literate.jl
# TODO notebooks
tutorial_dir = joinpath(docs_base, "src", "tutorials")
tutorial_output_dir = joinpath(docs_base, "src", "generated", "tutorials")
literate_config = Dict("credit" => false)
tutorial_pages = map(readdir(tutorial_dir)) do literate_source_file
    Literate.markdown(joinpath(tutorial_dir, literate_source_file), tutorial_output_dir;
                      config=literate_config)

    name = literate_source_file[1:end-3]
    path = joinpath("generated", "tutorials", literate_source_file[1:end-3] * ".md")
    name => path
end

makedocs(;
         modules=[Korg],
         repo=Documenter.Remotes.GitHub("ajwheeler", "Korg.jl"),
         sitename="Korg",
         pages=["Quickstart" => "index.md"
                "Install" => "install.md"
                "Tutorials" => [
                    "tutorials.md",
                    tutorial_pages...
                ]
                "Public Functions" => "API.md"
                "FAQ" => "FAQ.md"
                "Developer Documentation" => "devdocs.md"],

         # ideally we would do link checking, but adam doesn't like it turning his badge red
         linkcheck=false, # check if any external links are broken
         # there are useful if linkchecking gets turned back on
         linkcheck_timeout=30, # default (10s) is too short for slow academic sites
         linkcheck_ignore=["https://marcs.astro.uu.se/"],
         authors="Adam Wheeler, Matthew Abruzzo, Andrew Casey, and collaborators",
         format=Documenter.HTML(; assets=["assets/favicon.ico"], size_threshold=500000),
         # check that every function that has a docstring is in the docs
         checkdocs=:all,)

deploydocs(; repo="github.com/ajwheeler/Korg.jl.git", devbranch="main")
