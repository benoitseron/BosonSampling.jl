### Interferometers ###
"""
Supertype to any concrete interferomter type such as [`UserDefinedInterferometer`](@ref),
[`RandHaar`](@ref), [`Fourier`](@ref),...
"""
abstract type Interferometer end

"""
    LossyInterferometer <: Interferometer

Interferometers with inclusion of loss: the real `Interferometer` has dimension `m_real * m_real` while we model it by a `2m_real * 2m_real` one where the last `m_real` modes are environment modes containing the lost photons.
"""
abstract type LossyInterferometer <: Interferometer end

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
    transmission_amplitude::Real
    U::Matrix
    m::Int
    BeamSplitter(transmission_amplitude::Real) = new(transmission_amplitude, beam_splitter(transmission_amplitude), 2)
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
    angle::Real
    U::Matrix
    m::Int
    Rotation(angle::Real) = new(angle, rotation_matrix(angle),2)
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
    phase::Real
    U::Matrix
    m::Int
    PhaseShift(phase::Real) = new(phase, phase_shift(phase), 2)
end
