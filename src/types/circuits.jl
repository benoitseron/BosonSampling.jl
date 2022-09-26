abstract type Circuit <: Interferometer end
abstract type CircuitElement <: Interferometer end
abstract type LossyCircuitElement <: Interferometer end

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
    BeamSplitter(transmission_amplitude::Float64)

Creates a beam-splitter with tunable transmissivity.

    Fields:
        - transmission_amplitude::Float
        - U::Matrix{Complex}
        - m::Int
"""
struct BeamSplitter <: CircuitElement
    transmission_amplitude::Real
    U::Matrix
    m::Int
    BeamSplitter(transmission_amplitude::Real) = new(transmission_amplitude, beam_splitter(transmission_amplitude), 2)
end

LossParameters(::Type{BeamSplitter}) = IsLossless()

"""
    PhaseShift(phase)

Creates a phase-shifter with parameter `phase`.

    Fields:
        - phase::Float64
        - m::Int
        - U::Matrix{ComplexF64}
"""
struct PhaseShift <: CircuitElement
    phase::Real
    U::Matrix
    m::Int

    PhaseShift(phase) = new(phase, exp(1im * phase) * ones((1,1)), 1)

end

LossParameters(::Type{PhaseShift}) = IsLossless()

#
# """
#     Rotation(angle::Float64)
#
# Creates a Rotation matrix with tunable angle.
#
#     Fields:
#         - angle::Float64
#         - U::Matrix{ComplexF64}
#         - m::Int
# """
# struct Rotation <: Interferometer
#     angle::Real
#     U::Matrix
#     m::Int
#     Rotation(angle::Real) = new(angle, rotation_matrix(angle),2)
# end
#
# """
#     PhaseShift(phase::Float64)
#
# Creates a phase-shifter with parameter `phase`.
#
#     Fields:
#         - phase::Float64
#         - m::Int
#         - U::Matrix{ComplexF64}
# """
# struct PhaseShift <: Interferometer
#     phase::Real
#     U::Matrix
#     m::Int
#     PhaseShift(phase::Real) = new(phase, phase_shift(phase), 2)
# end


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
    RandomPhaseShifter <: CircuitElement

Applies a RandomPhase shift to a single optical line according to a distribution `d`. For a uniform phase shift for instance

    d = Uniform(0, 2pi)

"""
struct RandomPhaseShifter <: CircuitElement

    U::Matrix
    m::Int
    d::Distribution

    function RandomPhaseShifter(d::Distribution)

        new(exp(1im * rand(d)) * ones((1,1)), 1, d)

    end
end

LossParameters(::Type{RandomPhaseShifter}) = IsLossless()

"""
    is_compatible(target_modes_in::ModeList, circuit::Circuit)

Checks compatibility of `target_modes_in` and `circuit`.
"""
function is_compatible(target_modes_in::ModeList, circuit::Circuit)
    if circuit.m != target_modes_in.m
        if circuit.m == 2*target_modes_in.m
            error("use add_element_lossy! instead")
        else
            @show circuit.m
            @show target_modes_in.m
            error("target_modes_in.m")
        end
    end
    true
end

"""
    add_element!(circuit::Circuit, interf::Interferometer; target_modes::Vector{Int})
    add_element!(circuit::Circuit, interf::Interferometer; target_modes_in::Vector{Int}, target_modes_out::Vector{Int})

Adds the circuit element `interf` that will be applied on `target_modes` to the `circuit`.
Will automatically update the unitary representing the circuit.

If giving a single target modes, assumes that they are the same for out and in
"""
function add_element!(circuit::Circuit, interf::Interferometer, target_modes_in::ModeList, target_modes_out::ModeList = target_modes_in)

    @argcheck is_compatible(target_modes_in, target_modes_out)
    @argcheck is_compatible(target_modes_in, circuit)

    push!(circuit.circuit_elements, interf)
    circuit.U == nothing ? circuit.U = Matrix{ComplexF64}(I, circuit.m, circuit.m) : nothing

    if interf.m == circuit.m
        circuit.U = interf.U * circuit.U
    else
        u = Matrix{ComplexF64}(I, circuit.m, circuit.m)
        for i in 1:size(interf.U)[1]
            for j in 1:size(interf.U)[2]

                u[target_modes_in.modes[i], target_modes_out.modes[j]] = interf.U[i,j]

            end
        end

        # @show pretty_table(circuit.U)
        # @show pretty_table(u)

        circuit.U *= u

        # @show pretty_table(circuit.U)
    end

end

# # if giving a single target modes, assumes that they are the same for out and in
# function add_element!(circuit::Circuit, interf::Interferometer; target_modes::Vector{Int})
#     println("temporarily disabled")
#     # add_element!(circuit, interf; target_modes_in = target_modes, target_modes_out = target_modes)
#
# end
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

function add_element_lossy!(circuit::LossyCircuit, interf::Interferometer, target_modes_in::ModeList, target_modes_out::ModeList = target_modes_in)

    # @warn "health checks commented"

    if !(LossParameters(typeof(interf)) == IsLossy())
        # convert to a LossyInterferometer any lossless element just for size requirements and consistency

        println("converting to lossy")
        interf = to_lossy(interf)

    end

    if target_modes_in.m != circuit.m_real

        println("unexpected length")

        if target_modes_in.m == 2*circuit.m_real

            @warn "target_modes given with size 2*interf.m_real, discarding last m_real mode info and using the convention that mode i is lost into mode i+m_real"

        else

            @show circuit.m
            @show target_modes_in.m
            error("invalid size of target_modes_in.m")

    end

    end

    # @show target_modes_in
    # @show lossy_target_modes(target_modes_in)


     add_element!(circuit, interf, lossy_target_modes(target_modes_in), lossy_target_modes(target_modes_out))

end
#
# function add_element_lossy!(circuit::LossyCircuit, interf::Interferometer, target_modes_in::ModeOccupation, target_modes_out::ModeOccupation = target_modes_in)
#
#     target_modes_in = target_modes_in.state
#     target_modes_out = target_modes_out.state
#
#     @show target_modes_in
#
#     add_element_lossy!(circuit, interf, target_modes_in, target_modes_out)
#
# end
#
#
# function add_element_lossy!(circuit::LossyCircuit, interf::Interferometer, target_modes_in::ModeList, target_modes_out::ModeList = target_modes_in)
#
#     target_modes_in = convert(ModeOccupation, target_modes_in)
#     target_modes_out = convert(ModeOccupation, target_modes_out)
#
#     @show target_modes_in
#
#     add_element_lossy!(circuit, interf, target_modes_in, target_modes_out)
#
# end


Base.show(io::IO, interf::Interferometer) = begin
    print(io, "Interferometer :\n\n", "Type : ", typeof(interf), "\n", "m : ", interf.m, "\n", "U : ", "\n", interf.U)
end

Base.show(io::IO, interf::UserDefinedInterferometer) = begin
     print(io, "Interferometer :\n\n", "Type : ", typeof(interf), "\n", "m : ", interf.m, "\nUnitary : \n")
     pretty_table(io, interf.U)
end
