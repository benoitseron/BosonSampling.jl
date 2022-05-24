import Pkg
Pkg.add("Documenter")
using Documenter, BosonSampling
push!(LOAD_PATH, "./src")

makedocs(
    source = "./src/",
    sitename = "BosonSampling.jl",
    modules = [BosonSampling],
    authors = "Benoit Seron, Antoine Restivo",
    format = Documenter.HTML(prettyurls=false, sidebar_sitename=false),
    doctest = false,
    pages = [
        "About" => "about.md",
        "Tutorials" => Any[
                "installation" => "tutorial/installation.md",
                "basic usage" => "tutorial/basic_usage.md",
                "samplers" => "tutorial/boson_samplers.md",
                "bunching" => "tutorial/bunching.md",
                "certification" => "tutorial/certification.md",
                "optimization" => "tutorial/optimization.md",
                "circuits" => "tutorial/circuits.md",
                "permanent conjectures" => "tutorial/permanent_conjectures.md"
                        ],
        "API" => Any[
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
                "permanent conjectures" => "functions/permanent_conjectures.md",
                "tools" => "functions/proba_tools.md",
                "samplers" => "functions/samplers.md",
                "scattering" => "functions/scattering.md",
                "special_matrices" => "functions/special_matrices.md",
                "visualization" => "functions/visualize.md"]
                ]
    ]
)

deploydocs(
    repo = "github.com/BosonSampling.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"],
)
