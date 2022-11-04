# an interesting experiment to run

# use η_thermalization for the reflectivities
# we look at the number of photons in the last bin
# this would require pseudo-photon number resolution, but as you see, not up to
# many photons
# if we can recover the plot, it is one of the ways to show a difference between
# bosonic (giving a thermal distribution) versus distinguishable


"""
η_thermalization(n)

Defines the transmissivities required for the thermalization scheme.
"""
η_thermalization(n) = [1-((i-1)/i)^2 for i in 2:n]

"""
partition_thermalization(m)
Defines the last mode, single mode subset for thermalization. This corresponds to the first mode of the interferometer with spatial bins (to be checked).
"""
partition_thermalization(m) = begin

    s1 = Subset(ModeList(m,m))
    s2 = Subset(first_modes(m-1,m))
    Partition([s1,s2])

end


"""
    partition_thermalization_loss(m)

Defines the last mode, single mode subset for thermalization. This corresponds to the first mode of the interferometer with spatial bins (to be checked). Loss modes included in the second subset
"""
partition_thermalization_loss(m) = begin

    s1 = Subset(ModeList(m,2m))
    s2 = Subset(first_modes(m-1,2m)+last_modes(m,2m))


    Partition([s1,s2])

end

"""
    partition_thermalization_pnr(m)

Defines the modes corresponding to the pseudo number resolution as subsets for thermalization. This corresponds to the first mode of the interferometer with spatial bins. Loss modes included in the last subset.
"""
partition_thermalization_pnr(m) = begin

    subsets = [Subset(ModeList(i,2m)) for i in n:m]
    #push!(subsets, Subset(last_modes(m,2m))) # loss subset

    Partition(subsets)

end

"""
    η_pnr(steps)

Gives chosen refectivities for pseudo photon number resolution at the end of the loop with `steps`.
"""
η_pnr(steps) = [i/(steps + 1) for i in 1:steps]

"""

    tvd_reflectivities(η)

Computes the total variation distance for B,D inputs for a one loop setup with reflectivities η.
"""
function tvd_reflectivities(η)

    if all(between_one_and_zero.(η))

        params = LoopSamplingParameters(n=n, m=m, η = η, η_loss_bs = η_loss_bs, η_loss_lines = η_loss_lines, ϕ = ϕ)

        psp_b = convert(PartitionSamplingParameters, params)
        psp_d = convert(PartitionSamplingParameters, params) # a way to act as copy

        psp_b.mode_occ = equilibrated_input(sparsity, m)
        psp_d.mode_occ = equilibrated_input(sparsity, m)

        part = equilibrated_partition(m, n_subsets)

        psp_b.part = part
        psp_d.part = part

        psp_b.T = OneParameterInterpolation
        psp_b.x = x
        set_parameters!(psp_b)

        psp_d.T = Distinguishable
        set_parameters!(psp_d)

        compute_probability!(psp_b)
        compute_probability!(psp_d)

        pdf_bos = psp_b.ev.proba_params.probability.proba
        pdf_dist = psp_d.ev.proba_params.probability.proba

        @show tvd(pdf_bos,pdf_dist)
        #display(plot(η))
        return -tvd(pdf_bos,pdf_dist)
    else
        println("invalid reflectivity")
        return 100
    end

end
