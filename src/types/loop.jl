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
            if lossy && η_loss_lines != nothing
                interf = LossyLine(η_loss_lines[mode])
            else
                interf = RandomPhaseShifter(0.)
            end
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

    @unpack n, m, i, η, η_loss_bs, η_loss_lines, d, ϕ, p_dark, p_no_count = params

    build_loop(m, η, η_loss_bs, η_loss_lines, ϕ)

end


"""
    build_loop!(params::LoopSamplingParameters)
    build_loop!(data::OneLoopData)

Updates the `interferometer` field of a `LoopSamplingParameters`.
"""
function build_loop!(params::LoopSamplingParameters)
    params.interferometer = build_loop(params::LoopSamplingParameters)
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


