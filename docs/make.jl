using Documenter
using BosonSampling

makedocs(
    sitename = "BosonSampling",
    format = Documenter.HTML(),
    modules = [BosonSampling]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
