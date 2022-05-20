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
        "Library" => Any[
            "Types" => Any[
                "inputs" => "types/input.md",
                "events" => "types/events.md",
                "interferomters" => "types/interferometers.md",
                "measurements" => "types/measurements.md",
                "partitions" => "types/partitions.md",
                "type utilities" => "types/type_functions.md"],
            "Functions" => Any[
                "certification" => "functions/bayesian.md",
                "bunching" => "functions/bunching.md",
                "distributions" => "functions/distributions.md",
                "partitions" => "functions/partitions.md",
                "tools" => "functions/proba_tools.md",
                "samplers" => "functions/samplers.md",
                "scattering" => "functions/scattering.md",
                "special_matrices" => "functions/special_matrices.md"]
                ]
    ]
)

deploydocs(
    repo = "https://github.com/benoitseron/BosonSampling.jl.git"
)
