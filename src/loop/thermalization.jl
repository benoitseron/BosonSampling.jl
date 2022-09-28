include("packages_loop.jl")

# looking at the special interferometer of the work of antoine, with reflectivities set as below
# also tried with the random phases to see if it is affected

    n = 6
    m = n

    """
    η_thermalization(n)

    Defines the transmissivities required for the thermalization scheme.
    """
    η_thermalization(n) = [(i-1)/i for i in 2:n]


    """
    partition_thermalization(m)

    Defines the first mode, single mode subset for thermalization.
    """
    partition_thermalization(m) = begin

    s1 = Subset(ModeList(1,m))
    s2 = Subset(last_modes(m-1,m))
    Partition([s1,s2])

    end

    params = LoopSamplingParameters(n=n, η = η_thermalization(n), η_loss_bs = nothing, η_loss_lines = nothing, ϕ = nothing)

    @show η_thermalization(n)

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

    pdf_dist = pd.proba
    pdf_bos = pb.proba

    n_in = [i for i in 0:n]

    plot(n_in, pdf_bos, label = "B")
    plot!(n_in, pdf_dist, label = "D")

pb
