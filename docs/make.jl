push!(LOAD_PATH, "./src")

using Documenter, BosonSampling

About = "About" => "about.md"

Types = "Types" => "types.md"

PAGES = [About,Types]

makedocs(
    source = "./src/",
    sitename = "BosonSampling.jl",
    modules = [BosonSampling],
    authors = "Benoit Seron, Antoine Restivo",
    format = Documenter.HTML(),
    pages = PAGES
)

deploydocs(
    repo = "https://github.com/benoitseron/BosonSampling.jl.git"
)
