include("packages_loop.jl")

# looking at the special interferometer of the work of antoine, with reflectivities set as below
# also tried with the random phases to see if it is affected


"""
η_thermalization(n)

Defines the transmissivities required for the thermalization scheme.
"""
η_thermalization(n) = [(i-1)/i for i in 2:n]


"""
partition_thermalization(m)

Defines the last mode, single mode subset for thermalization. This corresponds to the first mode of the interferometer with spatial bins (to be checked).
"""
partition_thermalization(m) = begin

s1 = Subset(ModeList(m,m))
s2 = Subset(first_modes(m-1,m))
Partition([s1,s2])

end

@with_kw mutable struct PartitionSamplingParameters
    interf::Interferometer = RandHaar(m)
    T1::Type{T} where {T<:InputType} = Bosonic
    T2::Type{T} where {T<:InputType} = Distinguishable
    mode_occ_1::ModeOccupation = first_modes(n,m)
    mode_occ_2::ModeOccupation = first_modes(n,m)

    i1::Input = Input{T1}(mode_occ_1)
    i2::Input = Input{T1}(mode_occ_2)

    n_subsets::Int = 2
    part::Partition = equilibrated_partition(m, n_subsets)

    o::OutputMeasurementType = PartitionCountsAll(part)
    evb::Event = Event(ib,o,interf)
    evd::Event = Event(id,o,interf)

    ############## finish this below:

    pb = compute_probability!(evb)
    pd = compute_probability!(evd)

    pdf_dist::Union{Nothing, Vector{Float}} = pd.proba
    pdf_bos::Union{Nothing, Vector{Float}} = pb.proba

end

PartitionCountsAll<:OutputMeasurementType

@with_kw mutable struct PartitionLoopSamplingParameters
    params::LoopSamplingParameters

    @unpack n, m, input_type, i, η, η_loss_bs, η_loss_lines, d, ϕ, p_dark, p_no_count = params

    interf::Interferometer = build_loop(params)

    ####### partition parameters...

    part::Partition = equilibrated_partition(m, n_subsets)

end

### no phase ###

    n = 10
    m = n

    d = Uniform(0,2pi)
    ϕ = nothing # rand(d,m)


    params = LoopSamplingParameters(n=n, η = η_thermalization(n), η_loss_bs = nothing, η_loss_lines = nothing, ϕ = ϕ)

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

    pdf_dist = pd.proba
    pdf_bos = pb.proba

    n_in = [i for i in 0:n]

    plot(n_in, pdf_bos, label = "B", xticks = n_in)
    plot!(n_in, pdf_dist, label = "D")
    ylims!((0,1))
