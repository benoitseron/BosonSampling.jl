using Documenter
using BosonSampling

About = "Introduction" => "intro.md"

# GettingStarted = "gettingstarted.md"

Functions = "Functions" => "functions.md"

#

# Examples = "Examples" => [

#         "examples/flux.md"

#     ]

#

# License = "License" => "license.md"

PAGES = [About,  Functions]



makedocs(

    sitename = "BosonSampling.jl",

    modules = [BosonSampling],

    authors = "Benoit Seron, Antoine Restivo",

    format = Documenter.HTML(),

    pages = PAGES

)

# deploydocs(repo = "github.com/Evizero/Augmentor.jl.git")

# Documenter can also automatically deploy documentation to gh-pages.

# See "Hosting Documentation" and deploydocs() in the Documenter manual

# for more information.

deploydocs(

    repo = "github.com/benoitseron/BosonSampling.jl.git"

)
