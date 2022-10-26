include("packages_loop.jl")


###### what do we want to do ######

# find the parameters of reflectivities that allow for the best TVD, or an interesting behaviour
# fix the thermalization example to have this interesting case
# be able to optimize on reflectivity parameters
# for this I would like to be able to encapsulate everything in PartitionSamplingParameters etc


###### development ######



params = LoopSamplingParameters(m=10)



"""
    get_partition_sampling_parameters(params::LoopSamplingParameters, T2::Type{T} where {T<:InputType} = Distinguishable)

Unpacks the `params` to obtain a `PartitionSamplingParameters` compatible with as interferometer the circuit induced by `params`.
"""
function get_partition_sampling_parameters(params::LoopSamplingParameters)

    @unpack n, m, input_type, i, η, η_loss_bs, η_loss_lines, d, ϕ, p_dark, p_no_count = params

    PartitionSamplingParameters(n=n, m=m, T= get_parametric_type(params.i)[1], )

end


psp = convert(PartitionSamplingParameters, params)

compute_probability!(psp)

###### to threshold ######


mo = ModeOccupation([2,1,0])

ModeOccupation([(mode >= 1 ? 1 : 0) for mode in mo.state])

part = partition_thermalization_pnr(m)



[(length(subset)) for subset in part.subsets] == ones(length(part.subsets))

to_threshold(psp_b.ev.proba_params.probability.counts[10])

length(part.subsets[1])


# change a MultipleCounts to threshold

mc = psp_b.ev.proba_params.probability

typeof(mc)

mc.counts[1]



BosonSampling.to_threshold(mc::MultipleCounts)
