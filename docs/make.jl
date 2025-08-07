using Documenter, Literate, Suppressor, Korg # Korg is Pkg.dev'ed in the docs environment

using Pkg
Pkg.status()

# Check if we're in docs/ directory or root directory and set this path appropriately
# Locally, I often build from docs, but the CI builds from root
docs_base = basename(pwd()) == "docs" ? "." : "./docs"

# Use the README.md as index.md in the docs.
readme_path = joinpath(docs_base, "..", "README.md")
target_path = joinpath(docs_base, "src", "index.md")
cp(readme_path, target_path; force=true)

# paths for tutorial generation
tutorial_dir = joinpath(docs_base, "src", "Tutorials")
tutorial_output_dir = joinpath(docs_base, "src", "generated", "tutorials")
# delete everything in tutorial_output_dir
rm(tutorial_output_dir; force=true, recursive=true)
mkpath(tutorial_output_dir)

# process the tutorials into markdown and jupyter notebooks with Literate.jl
literate_config = Dict("credit" => false)
tutorial_pages = map(readdir(tutorial_dir)) do literate_source_file
    inpath = joinpath(tutorial_dir, literate_source_file)

    # generate markdown and notebook files
    markdown_path = Literate.markdown(inpath, tutorial_output_dir; config=literate_config)
    Literate.notebook(inpath, tutorial_output_dir; config=literate_config)

    # Make path relative to docs/src/. This is what Documenter.makedocs wants.
    rel_markdown_path = replace(markdown_path, r"^.*?docs/src/" => "")

    # page name => path
    literate_source_file[1:end-3] => rel_markdown_path
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

deploydocs(; repo="github.com/ajwheeler/Korg.jl.git",
           devbranch="main",
           push_preview=true) # https://ajwheeler.github.io/Korg.jl/previews/PR##
