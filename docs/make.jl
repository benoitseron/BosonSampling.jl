push!(LOAD_PATH, "path/to/BosonSampling.jl/src")

using Documenter, BosonSampling

About = "About" => " About.md"

GettingStarted = "gettingstarted.md"

Types = "Types" => "index.md"

PAGES = [About,Types]

makedocs(
    source = "/Path/to/BosonSampling.jl/docs/src/",
    sitename = "BosonSampling.jl",
    modules = [BosonSampling],
    authors = "Benoit Seron, Antoine Restivo",
    format = Documenter.HTML(),
    pages = PAGES
)

deploydocs(
    repo = "https://github.com/benoitseron/BosonSampling.jl.git"
)
