abstract type Loop <: Circuit end

"""
    LosslessCircuit(m::Int)

Creates an empty circuit with `m` input modes. The unitary representing the circuit
is accessed via the field `.U`.

    Fields:
        - m::Int
        - circuit_elements::Vector{Interferometer}
        - U::Union{Matrix{ComplexF64}, Nothing}
"""
mutable struct LosslessLoop <: Loop

    m::Int
    circuit_elements::Vector{Interferometer}
    U::Union{Matrix, Nothing}

    function LosslessLoop(m::Int)
        new(m, [], nothing)
    end

    function LosslessLoop(circuit_elements::Vector{Interferometer})

        new(length(circuit_elements) + 1, circuit_elements, nothing)
    end

end

LossParameters(::Type{LosslessLoop}) = IsLossless()

"""
    build_loop(m::Int, η::Union{T, Vector{T}}, η_loss_bs::Union{Nothing, T, Vector{T}} = nothing, η_loss_lines::Union{Nothing, T, Vector{T}} = nothing, ϕ::Union{Nothing, T, Vector{T}} = nothing) where {T<:Real}

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
            interf = LossyLineWithRandomPhase(η_loss_lines[mode], ϕ[mode])
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

            interf = LossyBeamSplitter(η[mode], η_loss_bs[mode])
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
