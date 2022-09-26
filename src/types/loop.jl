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

function build_loop(m::Int, η::Union{T, Vector{T}}, η_loss::Union{Nothing, T, Vector{T}} = nothing; lossy = false) where {T<:Real}

    if length(η) != m-1
        if length(η) == 1
            @warn "only a single transmissivity given - acts as a single beam splitter and not an array of beam splitters with a similar transmissivity"
        else
            error("invalid η, needs to be of size $(m-1)")
        end
    end

    if lossy == true
        error("not implemented")
    else
        if η_loss != nothing
            @warn "loss given but lossless loop"
        end

        circuit = LosslessCircuit(m)

        for mode in 1:m-1

            interf = BeamSplitter(η[mode])
            target_modes_in = [mode, mode+1]
            target_modes_out = [mode, mode+1]
            add_element!(circuit, interf, target_modes_in, target_modes_out)

        end

        return circuit.U
    end

end
