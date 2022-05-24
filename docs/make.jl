import Pkg
Pkg.add("Documenter")
using Documenter, BosonSampling
push!(LOAD_PATH, "./src")

DocMeta.setdocmeta!(BosonSampling, :DocTestSetup, :(using MyPackage); recursive=true)

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
                "Installation" => "tutorial/installation.md",
                "Basic usage" => "tutorial/basic_usage.md",
                "Samplers" => "tutorial/boson_samplers.md",
                "Bunching" => "tutorial/bunching.md",
                "Full distribution" => "tutorial/compute_distr.md",
                "Certification" => "tutorial/certification.md",
                "Optimization" => "tutorial/optimization.md",
                "Circuits" => "tutorial/circuits.md",
                "Permanent conjectures" => "tutorial/permanent_conjectures.md"
                        ],
        "API" => Any[
            "Types" => Any[
                "Inputs" => "types/input.md",
                "Events" => "types/events.md",
                "Interferomters" => "types/interferometers.md",
                "Measurements" => "types/measurements.md",
                "Partitions" => "types/partitions.md",
                "Type utilities" => "types/type_functions.md"],
            "Functions" => Any[
                "Certification" => "functions/bayesian.md",
                "Bunching" => "functions/bunching.md",
                "Distributions" => "functions/distributions.md",
                "Partitions" => "functions/partitions.md",
                "Permanent conjectures" => "functions/permanent_conjectures.md",
                "Tools" => "functions/proba_tools.md",
                "Samplers" => "functions/samplers.md",
                "Scattering" => "functions/scattering.md",
                "Special_matrices" => "functions/special_matrices.md",
                "Visualization" => "functions/visualize.md"]
                ]
    ]
)

deploydocs(
     repo = "github.com/BosonSampling.jl.git",
     target = "build",
     branch = "gh-pages",
     versions = ["stable" => "v^", "v#.#"],
 )
