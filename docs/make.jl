push!(LOAD_PATH, "./src")

using Documenter, BosonSampling

About = "About" => "about.md"

Types = "Types" => "types.md"

Functions = "Functions" => "functions.md"

PAGES = [About, Types, Functions]

makedocs(
    source = "./src/",
    sitename = "BosonSampling.jl",
    modules = [BosonSampling],
    authors = "Benoit Seron, Antoine Restivo",
    format = Documenter.HTML(prettyurls = false),
    pages = PAGES,
    Private = false
)

deploydocs(
    repo = "https://github.com/benoitseron/BosonSampling.jl.git"
)
