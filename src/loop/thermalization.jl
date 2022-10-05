include("packages_loop.jl")

# looking at the special interferometer of the work of antoine, with reflectivities set as below
# also tried with the random phases to see if it is affected


"""
η_thermalization(n)

Defines the transmissivities required for the thermalization scheme.
"""
η_thermalization(n) = [(i-1)/i for i in 2:n]






PartitionSamplingParameters(n = 10, m = 10)


PartitionCountsAll<:OutputMeasurementType

@with_kw mutable struct PartitionLoopSamplingParameters
    params::LoopSamplingParameters

    @unpack n, m, input_type, i, η, η_loss_bs, η_loss_lines, d, ϕ, p_dark, p_no_count = params

    interf::Interferometer = build_loop(params)

    ####### partition parameters...

    part::Partition = equilibrated_partition(m, n_subsets)

end


### no phase ###

begin

    n = 10
    m = n


    params = LoopSamplingParameters(n=n, η = η_thermalization(n), η_loss_bs = nothing, η_loss_lines = nothing) #, ϕ = nothing

    @show params

    @show η_thermalization(n)
    @show partition_thermalization(m)

    ib = Input{Bosonic}(first_modes(n,m))
    id = Input{Distinguishable}(first_modes(n,m))

    interf = build_loop(params)
    pretty_table(interf.U)
    part = partition_thermalization(m)
    o = PartitionCountsAll(part)
    evb = Event(ib,o,interf)
    evd = Event(id,o,interf)

    pb = compute_probability!(evb)
    pd = compute_probability!(evd)

    @show pdf_dist

    pdf_dist = pd.proba
    pdf_bos = pb.proba

    n_in = [i for i in 0:n]

    plot(n_in, pdf_bos, label = "B", xticks = n_in)
    plot!(n_in, pdf_dist, label = "D")
    ylims!((0,1))

end

typeof(pb)
