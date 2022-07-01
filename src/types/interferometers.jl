### Interferometers ###
"""
Supertype to any concrete interferomter type such as [`UserDefinedInterferometer`](@ref),
[`RandHaar`](@ref), [`Fourier`](@ref),...
"""
abstract type Interferometer end

"""
    UserDefinedInterferometer(U::Matrix)

Creates an instance of [`Interferometer`](@ref) from a given unitary matrix `U`.

    Fields:
        - m::Int
        - U::Matrix
"""
struct UserDefinedInterferometer <: Interferometer

    #actively checks unitarity, inefficient if outputing many matrices that are known to be unitary
    m::Int
    U::Matrix
    UserDefinedInterferometer(U) = is_unitary(U) ? new(size(U,1), U) : error("input matrix is non unitary")
end

"""
    RandHaar(m::Int)

Creates an instance of [`Interferometer`](@ref) from a Haar distributed unitary matrix of dimension `m`.

    Fields:
        - m::Int
        - U::Matrix
"""
struct RandHaar <: Interferometer
    m::Int
    U::Matrix{ComplexF64}
    RandHaar(m) = new(m,rand_haar(m))
end

"""
    Fourier(m::Int)

Creates a Fourier [`Interferometer`](@ref) of dimension `m`.

    Fields:
        - m::Int
        - U::Matrix
"""
struct Fourier <: Interferometer
    m::Int
    U::Matrix{ComplexF64}
    Fourier(m::Int) = new(m,fourier_matrix(m))
end

"""
    Hadamard(m::Int)

Creates a Hadamard [`Interferometer`](@ref) of dimension `m`.

    Fields:
        - m::Int
        - U::Matrix
"""
struct Hadamard <: Interferometer
    m::Int
    U::Matrix
    Hadamard(m::Int) = new(m,hadamard_matrix(m))
end

"""
    BeamSplitter(transmission_amplitude::Float64)

Creates a beam-splitter with tunable transmissivity.

    Fields:
        - transmission_amplitude::Float64
        - U::Matrix{ComplexF64}
        - m::Int
"""
struct BeamSplitter <: Interferometer
    transmission_amplitude::Float64
    U::Matrix
    m::Int
    BeamSplitter(transmission_amplitude::Float64) = new(transmission_amplitude, beam_splitter(transmission_amplitude), 2)
end

"""
    Rotation(angle::Float64)

Creates a Rotation matrix with tunable angle.

    Fields:
        - angle::Float64
        - U::Matrix{ComplexF64}
        - m::Int
"""
struct Rotation <: Interferometer
    angle::Float64
    U::Matrix
    m::Int
    Rotation(angle) = new(angle, rotation_matrix(angle),2)
end

"""
    PhaseShift(phase::Float64)

Creates a phase-shifter with parameter `phase`.

    Fields:
        - phase::Float64
        - m::Int
        - U::Matrix{ComplexF64}
"""
struct PhaseShift <: Interferometer
    phase::Float64
    U::Matrix
    m::Int
    PhaseShift(phase::Float64) = new(phase, phase_shift(phase), 2)
end

"""
    Circuit(m::Int)

Creates an empty circuit with `m` input modes. The unitary representing the circuit
is accessed via the field `.U`.

    Fields:
        - m::Int
        - circuit_elements::Vector{Interferometer}
        - U::Union{Matrix{ComplexF64}, Nothing}
"""
mutable struct Circuit <: Interferometer

    m::Int
    circuit_elements::Vector{Interferometer}
    U::Union{Matrix{ComplexF64}, Nothing}

    function Circuit(m::Int)
        new(m, [], nothing)
    end

end

"""
    add_element!(circuit::Circuit, interf::Interferometer, target_modes::Vector{Int})

Adds the circuit element `interf` that will be applied on `target_modes` to the `circuit`.
Will automatically update the unitary representing the circuit.      
"""
function add_element!(circuit::Circuit, interf::Interferometer, target_modes::Vector{Int})

    @argcheck interf.m == length(target_modes)

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
        circuit.U = u
    end

end

Base.show(io::IO, interf::Interferometer) = print(io, "Interferometer :\n\n", "Type : ", typeof(interf), "\n", "m : ", interf.m)

Base.show(io::IO, interf::UserDefinedInterferometer) = begin
     print(io, "Interferometer :\n\n", "Type : ", typeof(interf), "\n", "m : ", interf.m, "\nUnitary : \n")
     pretty_table(io, interf.U)
end
