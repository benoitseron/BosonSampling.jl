### Interferometers ###
"""
Supertype to any concrete interferomter type such as [`UserDefinedInterferometer`](@ref),
[`RandHaar`](@ref), [`Fourier`](@ref),...
"""
abstract type Interferometer end

"""
    IsLossy{T}
    IsLossless{T}

Trait to refer to quantities with inclusion of loss: the real `Interferometer` has dimension `m_real * m_real` while we model it by a `2m_real * 2m_real` one where the last `m_real` modes are environment modes containing the lost photons.
"""
abstract type LossParameters end

struct IsLossy <: LossParameters end
struct IsLossless <: LossParameters end

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

LossParameters(::Type{UserDefinedInterferometer}) = IsLossless()

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

LossParameters(::Type{RandHaar}) = IsLossless()

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

LossParameters(::Type{Fourier}) = IsLossless()

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

LossParameters(::Type{Hadamard}) = IsLossless()
