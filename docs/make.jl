push!(LOAD_PATH, "/Users/antoinerestivo/Desktop/BosonSampling.jl-optimized_types/src/")

using Documenter, BosonSampling

About = "About" => " About.md"

GettingStarted = "gettingstarted.md"

Types = "Types" => "index.md"

PAGES = [About,Types]

makedocs(
    source = "/Users/antoinerestivo/Desktop/BosonSampling.jl-optimized_types/docs/src/",
    sitename = "BosonSampling.jl",
    modules = [BosonSampling],
    authors = "Benoit Seron, Antoine Restivo",
    format = Documenter.HTML(),
    pages = PAGES
)

deploydocs(
    repo = "https://github.com/benoitseron/BosonSampling.jl.git"
)
