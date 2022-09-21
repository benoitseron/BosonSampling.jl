abstract type Circuit <: Interferometer end
abstract type LossyCircuitElement <: Circuit end

"""
    LosslessCircuit(m::Int)

Creates an empty circuit with `m` input modes. The unitary representing the circuit
is accessed via the field `.U`.

    Fields:
        - m::Int
        - circuit_elements::Vector{Interferometer}
        - U::Union{Matrix{ComplexF64}, Nothing}
"""
mutable struct LosslessCircuit <: Circuit

    m::Int
    circuit_elements::Vector{Interferometer}
    U::Union{Matrix, Nothing}

    function LosslessCircuit(m::Int)
        new(m, [], nothing)
    end

end

LossParameters(::Type{LosslessCircuit}) = IsLossless()

"""
    LossyCircuit(m_real::Int)

Lossy `Interferometer` constructed from `circuit_elements`.
"""
mutable struct LossyCircuit <: Circuit

    m_real::Int
    m::Int
    circuit_elements::Vector{Interferometer}
    U_physical::Union{Matrix{Complex}, Nothing} # physical matrix
    U::Union{Matrix{Complex}, Nothing} # virtual 2m*2m interferometer

    function LossyCircuit(m_real::Int)
        new(m_real, 2*m_real, [], nothing, nothing)
    end

end

LossParameters(::Type{LossyCircuit}) = IsLossy()

"""
    LossyBeamSplitter(transmission_amplitude, η_loss)

Creates a beam-splitter with tunable transmissivity and loss. Uniform model of loss: each input line i1,i2 has a beam splitter in front with transmission amplitude of `transmission_amplitude_loss` into an environment mode.
"""
struct LossyBeamSplitter <: LossyCircuitElement
    transmission_amplitude::Real
    η_loss::Real
    U::Matrix
    m::Int
    m_real::Int
    function LossyBeamSplitter(transmission_amplitude::Real, η_loss::Real)
        @argcheck between_one_and_zero(transmission_amplitude)
        @argcheck between_one_and_zero(η_loss)

        new(transmission_amplitude, η_loss, virtual_interferometer_uniform_loss(beam_splitter(transmission_amplitude),η_loss), 4,2)
    end
end

LossParameters(::Type{LossyBeamSplitter}) = IsLossy()

"""
    LossyLine(η_loss)

Optical line with some loss: represented by a `BeamSplitter` with a transmission amplitude of `transmission_amplitude_loss` into an environment mode.
"""
struct LossyLine <: LossyCircuitElement

    η_loss::Real
    U::Matrix
    m::Int
    m_real::Int

    function LossyLine(η_loss::Real)

        @argcheck between_one_and_zero(η_loss)

        new(η_loss, virtual_interferometer_uniform_loss(ones((1,1)), η_loss), 2,1)
    end
end

LossParameters(::Type{LossyLine}) = IsLossy()

"""
    add_element!(circuit::Circuit, interf::Interferometer; target_modes::Vector{Int})
    add_element!(;circuit::Circuit, interf::Interferometer, target_modes::Vector{Int})

Adds the circuit element `interf` that will be applied on `target_modes` to the `circuit`.
Will automatically update the unitary representing the circuit.
"""
function add_element!(circuit::Circuit, interf::Interferometer; target_modes::Vector{Int})

    if interf.m != length(target_modes)
        if interf.m == 2*length(target_modes)
            error("use add_element_lossy! instead")
        else
            error("invalid length(target_modes)")
        end
    end

    push!(circuit.circuit_elements, interf)
    circuit.U == nothing ? circuit.U = Matrix{ComplexF64}(I, circuit.m, circuit.m) : nothing

    if interf.m == circuit.m
        circuit.U = interf.U * circuit.U
    else
        u = Matrix{ComplexF64}(I, circuit.m, circuit.m)
        for i in 1:size(interf.U)[1]
            for j in 1:size(interf.U)[2]
                u[target_modes[i], target_modes[j]] = interf.U[i,j]
            end
        end

        @show circuit.U
        @show u

        circuit.U = u * circuit.U

        @show circuit.U
    end

end
#
# add_element!(circuit::Circuit, interf::Interferometer; target_modes::ModeOccupation) = add_element!(circuit=circuit, interf=interf, target_modes=target_modes.state)
#
# add_element!(circuit::Circuit, interf::Interferometer; target_modes::Vector{Int}) = add_element!(circuit=circuit, interf=interf, target_modes=target_modes)
#

# bug:

# WARNING: Method definition add_element!(BosonSampling.Circuit, BosonSampling.Interferometer) in module BosonSampling at /home/benoitseron/.julia/dev/BosonSampling/src/types/circuits.jl:73 overwritten at /home/benoitseron/.julia/dev/BosonSampling/src/types/circuits.jl:75.
#   ** incremental compilation may be fatally broken for this module **
#
# WARNING: Method definition add_element!##kw(Any, typeof(BosonSampling.add_element!), BosonSampling.Circuit, BosonSampling.Interferometer) in module BosonSampling at /home/benoitseron/.julia/dev/BosonSampling/src/types/circuits.jl:73 overwritten at /home/benoitseron/.julia/dev/BosonSampling/src/types/circuits.jl:75.
#   ** incremental compilation may be fatally broken for this module **

# add_element!(circuit::Circuit, interf::Interferometer; target_modes::ModeList) = begin
#
#     mo = convert(ModeOccupation, target_modes)
#     add_element!(circuit=circuit, interf=interf, target_modes=mo)
#
# end

function add_element_lossy!(circuit::LossyCircuit, interf::Interferometer; target_modes::Vector{Int})

    # @warn "health checks commented"

    if !(LossParameters(typeof(interf)) == IsLossy())
        # convert to a LossyInterferometer any lossless element just for size requirements and consistency

        println("converting to lossy")
        interf = to_lossy(interf)

    end

    if length(target_modes) != interf.m_real

        println("unexpected length")

        if length(target_modes) == 2*interf.m_real

            @warn "target_modes given with size 2*interf.m_real, discarding last m_real mode info and using the convention that mode i is lost into mode i+m_real"

        else

            error("invalid size of target_modes")

        end

    end

    # @show lossy_target_modes(target_modes)
     add_element!(circuit, interf, target_modes=lossy_target_modes(target_modes))

end


Base.show(io::IO, interf::Interferometer) = print(io, "Interferometer :\n\n", "Type : ", typeof(interf), "\n", "m : ", interf.m)

Base.show(io::IO, interf::UserDefinedInterferometer) = begin
     print(io, "Interferometer :\n\n", "Type : ", typeof(interf), "\n", "m : ", interf.m, "\nUnitary : \n")
     pretty_table(io, interf.U)
end
