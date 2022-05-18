push!(LOAD_PATH, "./src")

using Documenter, BosonSampling

makedocs(
    source = "./src/",
    sitename = "BosonSampling.jl",
    modules = [BosonSampling],
    authors = "Benoit Seron, Antoine Restivo",
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "About" => "about.md",
        "Types" => Any[
            "inputs" => "types/input.md",
            "events" => "types/events.md"],
        "Functions" => Any[
            "bunching" => "functions/bunching.md"
        ]
    ]
)

deploydocs(
    repo = "https://github.com/benoitseron/BosonSampling.jl.git"
)
