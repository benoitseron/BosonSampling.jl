include("packages_loop.jl")

"""
    η_pnr(steps)

Gives chosen refectivities for pseudo photon number resolution at the end of the loop with `steps`.
"""
η_pnr(steps) = [i/(steps + 1) for i in 1:steps]
