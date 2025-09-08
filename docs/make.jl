# the file is run by CI to build the docs.
# you can run it locally with `julia --project=docs/Project.toml make.jl` from the docs/ directory
# add --draft to skip executing code blocks and running tutorial notebooks for faster builds

using Documenter, Literate, Suppressor, ArgParse, Korg # Korg is Pkg.dev'ed in the docs environment

# Set up command line arguments
s = ArgParseSettings()
#! format: off
@add_arg_table! s begin
    "--draft"
    help = "Don't execute code blocks or run tutorial notebooks. Run makedocs in draft mode."
    action = :store_true
    default = false
end
#! format: on
parsed_args = parse_args(s)

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
if parsed_args["draft"]
    # don't change if --no-execute is not passed, because the defaults for md and nb are different
    literate_config["execute"] = false
end
tutorial_pages = map(readdir(tutorial_dir)) do literate_source_file
    inpath = joinpath(tutorial_dir, literate_source_file)

    # generate markdown and notebook files
    markdown_path = Literate.markdown(inpath, tutorial_output_dir; config=literate_config)
    @time Literate.notebook(inpath, tutorial_output_dir; config=literate_config)

    # Make path relative to docs/src/. This is what Documenter.makedocs wants.
    rel_markdown_path = replace(markdown_path, r"^.*?docs/src/" => "")

    # page name => path
    literate_source_file[1:end-3] => rel_markdown_path
end

makedocs(;
         modules=[Korg],
         doctest=!parsed_args["draft"],
         draft=parsed_args["draft"],
         repo=Documenter.Remotes.GitHub("ajwheeler", "Korg.jl"),
         sitename="Korg",
         pages=["Quickstart" => "index.md"
                "Install" => "install.md"
                "Upgrading to v1.0" => "Upgrading.md"
                "Guides" => [
                    tutorial_pages...,
                    "Other Tutorials" => "tutorials.md",
                    "Abundances.md",
                    "Wavelengths.md"
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
