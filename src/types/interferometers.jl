### Interferometers ###

abstract type Interferometer end

"""
    UserDefinedInterferometer(U::Matrix)

Creates an instance of `Interferometer` from a provided unitary matrix ``U``.
"""
struct UserDefinedInterferometer <: Interferometer

    #actively checks unitarity, inefficient if outputing many matrices that are known to be unitary
    m::Int
    U::Matrix
    UserDefinedInterferometer(U) = is_unitary(U) ? new(size(U,1), U) : error("input matrix is non unitary")
end

"""
    RandHaar(m::Int)

Creates an instance of `Interferometer` from a Haar distributed unitary matrix of dimension ``m``.

!!! tip "Fields"
    - `m::Int`: Matrix dimension,
    - `U::Matrix{ComplexF64}`: Matrix representation of the interferometer
"""
struct RandHaar <: Interferometer
    m::Int
    U::Matrix{ComplexF64}
    RandHaar(m) = new(m,rand_haar(m))
end

"""
    Fourier(m::Int)

Creates a Fourier `Interferometer` of dimension ``m``.

!!! tip "Fields"
    - `m::Int`: Matrix dimension,
    - `U::Matrix{ComplexF64}`: Matrix representation of the interferometer
"""
struct Fourier <: Interferometer
    m::Int
    U::Matrix{ComplexF64}
    Fourier(m::Int) = new(m,fourier_matrix(m))
end

"""
    Hadamard(m::Int)

Creates a Hadamard `Interferometer` of dimension ``m``.

!!! tip "Fields"
    - `m::Int`: Matrix dimension,
    - `U::Matrix{ComplexF64}`: Matrix representation of the interferometer
"""
struct Hadamard <: Interferometer
    m::Int
    U::Matrix
    Hadamard(m::Int) = new(m,hadamard_matrix(m))
end

"""
    BeamSplitter(transmission_amplitude::Float64)

Creates a beam-splitter with tunable transmissivity.

!!! tip "Fields"
    - `transmission_amplitude::Float64`,
    - `U::Matrix{ComplexF64}`: Matrix representation of the interferometer,
    - `m::Int`: Matrix dimension
"""
struct BeamSplitter <: Interferometer
    transmission_amplitude::Float64
    U::Matrix
    m::Int
    BeamSplitter(transmission_amplitude) = new(transmission_amplitude, beam_splitter(transmission_amplitude),2)
end

"""
    Rotation(angle::Float64)

Creates a Rotation matrix with tunable angle.

!!! tip "Fields"
    - `angle::Float64`,
    - `U::Matrix{ComplexF64}`: Matrix representation of the interferometer,
    - `m::Int`: Matrix dimension
"""
struct Rotation <: Interferometer
    angle::Float64
    U::Matrix
    m::Int
    Rotation(angle) = new(angle, rotation_matrix(angle),2)
end

"""
    PhaseShift(shifted_modes::Array, param_::Array)

Creates a phase-shifter that is applied on the modes precised by shifted_modes with phase shifts given in param_.

!!! tip "Fields"
    - `shifted_modes::Array`,
    - `param_::Array`,
    - `m::Int`: Matrix dimension,
    - `U::Matrix{ComplexF64}`: Matrix representation of the interferometer
"""
struct PhaseShift <: Interferometer
    shifted_modes::Array
    param_::Array
    m::Int
    U::Matrix
    PhaseShift(shifted_modes, param_) = new(shifted_modes, param_, length(shifted_modes), phase_shift(shifted_modes, param_))
end

Base.show(io::IO, interf::Interferometer) = print(io, "Interferometer :\n\n", "Type : ", typeof(interf), "\n", "m : ", interf.m)

Base.show(io::IO, interf::UserDefinedInterferometer) = begin
     print(io, "Interferometer :\n\n", "Type : ", typeof(interf), "\n", "m : ", interf.m, "\nUnitary : \n")
     pretty_table(io, interf.U)
end
