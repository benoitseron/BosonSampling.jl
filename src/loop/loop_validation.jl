include("packages_loop.jl")

###### fixed random unitary ######

# fix a random unitary with potential noise (as of now just rand phases)
# find the partition data and tvd
# repeat over niter
# reflectivities chosen at random

"""
    loop_partition_tvd(params::LoopSamplingParameters, n_subsets::Int = 2)

Outputs the TVD between the partition probabilities for `Bosonic`, `Distinguishable` over the Interferometer defined by `params` (builds a new loop at each call so if they contain randomness this is taken into account however if the randomness is called at the definition of params you need to rerun it before feeding it into the function).
"""
function loop_partition_tvd(params::LoopSamplingParameters, n_subsets::Int = 2)

    ib = Input{Bosonic}(first_modes(n,m))
    id = Input{Distinguishable}(first_modes(n,m))

    interf = build_loop(params)
    part = equilibrated_partition(m,n_subsets)
    o = PartitionCountsAll(part)
    evb = Event(ib,o,interf)
    evd = Event(id,o,interf)

    pb = compute_probability!(evb)
    pd = compute_probability!(evd)

    pdf_dist = pd.proba
    pdf_bos = pb.proba

    tvd(pdf_bos,pdf_dist)

end

begin
    n = 10
    m = n
    niter = 1000
    n_subsets = 2

    tvd_array = zeros(niter)

    for i in 1:niter

        params = LoopSamplingParameters(n=n, η_loss_bs = nothing, η_loss_lines = nothing)

        tvd_array[i] = loop_partition_tvd(params, n_subsets)
    end

    mean(tvd_array), sqrt.(var(tvd_array))

end
# end

params = LoopSamplingParameters(n = n, input_type = Distinguishable, η_loss_lines = 0.3 * ones(m))

params.η_loss_lines

build_loop(params).U

pretty_table(build_loop(params).U)

n = 2
m = n
get_sample_loop(LoopSamplingParameters(n = n, input_type = Distinguishable, η_loss_lines = 0.3 * ones(m)))
