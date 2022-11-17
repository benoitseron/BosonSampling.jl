
# """
#         PartitionSamplingParameters
#
# Type holding info on a numerical experiment simulating the probability distributions of photon counting in partitions, used for comparing two input types `T1`,`T2`.
#
# By default the interferometer is set at nothing - you need to use `set_interferometer!`.
#
#         Fields:
#                 n::Int
#                 m::Int
#                 interf::Interferometer = RandHaar(m)
#
#                 T1::Type{T} where {T<:InputType} = Bosonic
#                 T2::Type{T} where {T<:InputType} = Distinguishable
#                 mode_occ_1::ModeOccupation = first_modes(n,m)
#                 mode_occ_2::ModeOccupation = first_modes(n,m)
#
#                 i1::Input = Input{T1}(mode_occ_1)
#                 i2::Input = Input{T1}(mode_occ_2)
#
#                 n_subsets::Int = 2
#                 part::Partition = equilibrated_partition(m, n_subsets)
#
#                 o::OutputMeasurementType = PartitionCountsAll(part)
#                 ev1::Event = Event(i1,o,interf)
#                 ev2::Event = Event(i2,o,interf)
# """
# @with_kw mutable struct PartitionSamplingParameters
#
#     n::Int
#     m::Int = n
#     m_real::Int = m
#     interf::Union{Interferometer, Nothing} = nothing
#
#     # begin
#     #     if interf != nothing
#     #         @warn "use set_interferometer! to update the events as well as everything in a lossy form (if it is a lossy interferometer)"
#     #     end
#     # end
#
#     T1::Type{T} where {T<:InputType} = Bosonic
#     T2::Type{T} where {T<:InputType} = Distinguishable
#     mode_occ_1::ModeOccupation = first_modes(n,m)
#     mode_occ_2::ModeOccupation = first_modes(n,m)
#
#     i1::Input = Input{T1}(mode_occ_1)
#     i2::Input = Input{T2}(mode_occ_2)
#
#     n_subsets::Int = 2
#
#     part::Partition = equilibrated_partition(m, n_subsets)
#
#     o::OutputMeasurementType = PartitionCountsAll(part)
#     ev1::Union{Event, Nothing} = nothing
#     ev2::Union{Event, Nothing} = nothing
#
#
#
# end

@with_kw mutable struct SamplingParameters

    n::Int
    m::Int = n
    m_real::Int = m
    interf::Union{Interferometer, Nothing} = nothing

    T::Type{T} where {T<:InputType} = Bosonic
    mode_occ::ModeOccupation = first_modes(n,m)
    x::Union{Nothing, Real} = nothing

    i::Union{Input, Nothing} = nothing

    o::Union{OutputMeasurementType, Nothing} = nothing
    ev::Union{Event, Nothing} = nothing

end


@with_kw mutable struct PartitionSamplingParameters

    n::Int
    m::Int = n
    m_real::Int = m
    interf::Union{Interferometer, Nothing} = nothing

    # begin
    #     if interf != nothing
    #         @warn "use set_interferometer! to update the events as well as everything in a lossy form (if it is a lossy interferometer)"
    #     end
    # end

    T::Type{T} where {T<:InputType} = Bosonic
    mode_occ::ModeOccupation = first_modes(n,m)
    x::Union{Nothing, Real} = nothing

    i::Union{Input, Nothing} = nothing

    # if T == OneParameterInterpolation
    #     i = Input{T1}(mode_occ,x)
    # elseif T in [Bosonic, Distinguishable]
    #     i = Input{T1}(mode_occ)
    # else
    #     error("type not implemented")
    # end

    n_subsets::Int = 2

    part::Partition = equilibrated_partition(m, n_subsets)

    o::OutputMeasurementType = PartitionCountsAll(part)
    ev::Union{Event, Nothing} = nothing

end

"""
    set_input!(params::PartitionSamplingParameters)
    set_input!(params::SamplingParameters)

`PartitionSamplingParameters` is initially defined with a nothing input, this function acts as an outer constructor to make it compatible with peculiarities of the `@with_kw` used in the type definition.
"""
function set_input!(params::Union{PartitionSamplingParameters, SamplingParameters})

    if params.T in [Bosonic, Distinguishable]
        params.i = Input{params.T}(params.mode_occ)
    elseif params.T == OneParameterInterpolation
        params.i = Input{params.T}(params.mode_occ, params.x)
    else
        error("type not implemented")
    end

end

"""
    set_interferometer!(interf::Interferometer, params::PartitionSamplingParameters)
    function set_interferometer!(params::PartitionSamplingParameters)

Updates the interferometer in a PartitionSamplingParameters, including the definition of the events and upgrade to lossy if needed of other quantities.

This function acts as an outer constructor to make it compatible with peculiarities of the `@with_kw` used in the type definition.
"""
function set_interferometer!(interf::Interferometer, params::Union{PartitionSamplingParameters, SamplingParameters})

    if params.i == nothing
        set_input!(params)
    end

    params.interf = interf

    if LossParameters(typeof(interf)) == IsLossy()

        @argcheck interf.m == 2*params.m

        for field in [:mode_occ, :i,:part, :o]
            getfield(params, field) = to_lossy(getfield(params, field))
        end

        if typeof(params) == PartitionSamplingParameters
            params.n_subsets += 1
        end
    else
        @argcheck interf.m == params.m
    end

    params.ev =  Event(params.i,params.o,params.interf)

end

set_interferometer!(params::Union{PartitionSamplingParameters, SamplingParameters}) = set_interferometer!(params.interf, params)

function set_partition!(params::PartitionSamplingParameters)

    params.o = PartitionCountsAll(params.part)
    params.ev =  Event(params.i,params.o,params.interf)

end

function set_measurement!(o::OutputMeasurementType, params::SamplingParameters)


    if StateMeasurement(typeof(o)) == FockStateMeasurement()
        params.o = o
        params.ev =  Event(params.i,params.o,params.interf)
    else
        error("invalid measurement")

    end

end

set_measurement!(params::SamplingParameters) = set_measurement!(params.o, params)

function set_parameters!(params::Union{PartitionSamplingParameters, SamplingParameters})

    set_input!(params)

    if typeof(params) == PartitionSamplingParameters
        set_partition!(params)
    elseif typeof(params) == SamplingParameters
        set_measurement!(params)
    end

    set_interferometer!(params)

end


"""
    LoopSamplingParameters(...)

Container for sampling parameters with a LoopSampler. Parameters are set by default as defined, and you can change only the ones needed, for instance to sample `Distinguishable` particles instead, just do

    LoopSamplingParameters(input_type = Distinguishable)

and to change the number of photons with it

    LoopSamplingParameters(n = 10 ,input_type = Distinguishable)

To be used with [`get_sample_loop`](@ref).

By default it applies a random phase at each optical line.
"""
@with_kw mutable struct LoopSamplingParameters

    n::Int = 4
    m::Int = n
    x::Union{Real, Nothing} = nothing
    input_type::Type{T} where {T<:InputType} = Bosonic
    i::Input = begin
        if input_type in [Bosonic, Distinguishable]
            i =  Input{input_type}(first_modes(n,m))
        elseif input_type == OneParameterInterpolation
            if x == nothing
                error("x not given")
            else
                i = Input{input_type}(first_modes(n,m), x)
            end
        end
    end


    η::Union{T, Vector{T}}  where {T<:Real} = 1/sqrt(2) .* ones(m-1)
    η_loss_bs::Union{Nothing, T, Vector{T}}   where {T<:Real} = 1 .* ones(m-1)
    η_loss_lines::Union{Nothing, T, Vector{T}}   where {T<:Real} = 1 .* ones(m)
    d::Union{Nothing, Real, Distribution} = Uniform(0, 2pi)
    ϕ::Union{Nothing, T, Vector{T}} where {T<:Real} = rand(d, m)

    p_dark::Real = 0.0
    p_no_count::Real = 0.0

end


function Base.convert(::Type{PartitionSamplingParameters}, params::LoopSamplingParameters)

    @unpack n, m, input_type, i, η, η_loss_bs, η_loss_lines, d, ϕ, p_dark, p_no_count = params

    interf = build_loop(params)

    ps = PartitionSamplingParameters(n=n, m=m, T= get_parametric_type(i)[1], interf = interf, mode_occ = i.r, x = i.distinguishability_param,i=i)

    set_interferometer!(ps)

    ps

end
