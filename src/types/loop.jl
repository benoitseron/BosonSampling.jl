abstract type Loop <: Circuit end

"""
    LosslessLoop(m, η, ϕ)

Creates a LosslessLoop, see [`build_loop`](@ref) for the fields.
"""
mutable struct LosslessLoop<: Loop

    m::Int
    η::Union{T, Vector{T}}  where {T<:Real}
    ϕ::Union{Nothing, T, Vector{T}}   where {T<:Real}
    circuit::Circuit
    U::Union{Nothing, Matrix}

    function LosslessLoop(m, η, ϕ)
        circuit = build_loop(m, η, nothing, nothing, ϕ)
        new(m, η, ϕ, circuit, circuit.U)

    end
end

LossParameters(::Type{LosslessLoop}) = IsLossless()

"""
    LossyLoop(m::Int)

Creates a LossyLoop, see [`build_loop`](@ref) for the fields.
"""
mutable struct LossyLoop<: Loop

    m::Int
    η::Union{T, Vector{T}}  where {T<:Real}
    η_loss_bs::Union{Nothing, T, Vector{T}} where {T<:Real}
    η_loss_lines::Union{Nothing, T, Vector{T}} where {T<:Real}
    ϕ::Union{Nothing, T, Vector{T}}   where {T<:Real}
    circuit::Circuit
    U::Union{Nothing, Matrix}

    function LossyLoop(m, η, η_loss_bs, η_loss_lines, ϕ)
        circuit =  build_loop(m, η, η_loss_bs, η_loss_lines, ϕ)
        new(m, η, η_loss_bs, η_loss_lines, ϕ, circuit, circuit.U)

    end
end

LossParameters(::Type{LossyLoop}) = IsLossy()


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
    input_type::Type{T} where {T<:InputType} = Bosonic
    i::Input = Input{input_type}(first_modes(n,m))

    η::Union{T, Vector{T}}  where {T<:Real} = 1/sqrt(2) .* ones(m-1)
    η_loss_bs::Union{Nothing, T, Vector{T}}   where {T<:Real} = 1 .* ones(m-1)
    η_loss_lines::Union{Nothing, T, Vector{T}}   where {T<:Real} = 1 .* ones(m)
    d::Union{Nothing, Real, Distribution} = Uniform(0, 2pi)
    ϕ::Union{Nothing, T, Vector{T}}   where {T<:Real} = rand(d, m)

    p_dark::Real = 0.0
    p_no_count::Real = 0.0

end

"""
    build_loop(m::Int, η::Union{T, Vector{T}}, η_loss_bs::Union{Nothing, T, Vector{T}} = nothing, η_loss_lines::Union{Nothing, T, Vector{T}} = nothing, ϕ::Union{Nothing, T, Vector{T}} = nothing) where {T<:Real}

    build_loop(params::LoopSamplingParameters)

Outputs a `Circuit` corresponding to the experiment discussed with the Walther group.

A pulse of `n` photons in `m` temporal modes is sent through a variable `BeamSplitter` and a `LossyLineWithRandomPhase` as described in "Scalable Boson Sampling with Time-Bin Encoding Using a Loop-Based Architecture".

    Fields:
        - m::Int dimension
        - η::Union{T, Vector{T}} array of beam splitter transmissivities
        - η_loss_bs::Union{Nothing, T, Vector{T}} array of beam splitter transmissivities for accounting loss
        - η_loss_lines::Union{Nothing, T, Vector{T}} array of beam delay lines transmissivities for accounting loss
        - ϕ::Union{Nothing, T, Vector{T}} array of phases applied by the delay lines
"""
function build_loop(m::Int, η::Union{T, Vector{T}}, η_loss_bs::Union{Nothing, T, Vector{T}} = nothing, η_loss_lines::Union{Nothing, T, Vector{T}} = nothing, ϕ::Union{Nothing, T, Vector{T}} = nothing) where {T<:Real}

    for param in [η, η_loss_bs, η_loss_lines, ϕ]
        if param != nothing
            if length(param) != m-1
                if length(param) == 1
                    @warn "only a single $(Symbol(param)) given - acts as a single beam splitter and not an array of beam splitters with a similar transmissivity"
                else
                    if param in [ϕ, η_loss_lines]
                        if length(param) != m
                            error("invalid $(Symbol(param)), needs to be of size $(m)")
                        end
                    else
                        error("invalid $(Symbol(param)), needs to be of size $(m-1)")
                    end
                end
            end
        end
    end

    function add_line!(mode, lossy)

        if ϕ == nothing
            error("not implemented")
        else
            if η_loss_lines == nothing
                interf = RandomPhaseShifter(ϕ[mode])
            else
                interf = LossyLineWithRandomPhase(η_loss_lines[mode], ϕ[mode])
            end

            target_modes_in = ModeList([mode], circuit.m_real)
            target_modes_out = target_modes_in

            if lossy
                add_element_lossy!(circuit, interf, target_modes_in, target_modes_out)
            else
                add_element!(circuit, interf, target_modes_in, target_modes_out)
            end
        end
    end

    if η_loss_bs != nothing ||  η_loss_lines != nothing

        lossy = true

        circuit = LossyCircuit(m)

        for mode in 1:m-1

            add_line!(mode, lossy)

            if η_loss_bs != nothing
                interf = LossyBeamSplitter(η[mode], η_loss_bs[mode])
            else
                interf = LossyBeamSplitter(η[mode], 1.)
            end

            target_modes_in = ModeList([mode, mode+1], circuit.m_real)
            target_modes_out = target_modes_in
            add_element_lossy!(circuit, interf, target_modes_in, target_modes_out)

        end

        add_line!(m, lossy) # last pass has no interaction with a BS

        return circuit

    else

        lossy = false

        circuit = LosslessCircuit(m)

        for mode in 1:m-1

            add_line!(mode, lossy)
            interf = BeamSplitter(η[mode])#LossyBeamSplitter(η[mode], η_loss[mode])
            #target_modes_in = ModeList([mode, mode+1], circuit.m_real)
            #target_modes_out = ModeList([mode, mode+1], circuit.m_real)

            target_modes_in = ModeList([mode, mode+1], m)
            target_modes_out = target_modes_in
            add_element!(circuit, interf, target_modes_in, target_modes_out)

        end
        add_line!(m, lossy)

        return circuit
    end

end

function build_loop(params::LoopSamplingParameters)

    @unpack n, m, input_type, i, η, η_loss_bs, η_loss_lines, d, ϕ, p_dark, p_no_count = params

    build_loop(m, η, η_loss_bs, η_loss_lines, ϕ)

end



"""
    get_sample_loop(params::LoopSamplingParameters)

Obtains a sample for `LoopSamplingParameters` by reconstructing the circuit each time (as is needed for adding a random phase).
"""
function get_sample_loop(params::LoopSamplingParameters)

    @unpack n, m, input_type, i, η, η_loss_bs, η_loss_lines, d, ϕ, p_dark, p_no_count = params

    circuit = LossyLoop(m, η, η_loss_bs, η_loss_lines, ϕ).circuit

    o = RealisticDetectorsFockSample(p_dark, p_no_count)
    ev = Event(i,o, circuit)

    BosonSampling.sample!(ev)

    ev.output_measurement.s

end
