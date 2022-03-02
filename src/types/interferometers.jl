### Interferometers ###

abstract type Interferometer end

struct UserDefinedInterferometer <: Interferometer

    #actively checks unitarity, inefficient if outputing many matrices that are known to be unitary
    m::Int
    U::Matrix
    UserDefinedInterferometer(U::Matrix) = is_unitary(U) ? new(m, U) : error("input matrix is non unitary")
end

struct RandHaar <: Interferometer
    m::Int
    U::Matrix{ComplexF64}
    RandHaar(m) = new(m,rand_haar(m))
end

struct Fourier <: Interferometer
    m::Int
    U::Matrix{ComplexF64}
    Fourier(m::Int) = new(m,fourier_matrix(m))
end


struct BeamSplitter <: Interferometer
    transmission_amplitude::Float64
    U::Matrix
    BeamSplitter(transmission_amplitude) = new(transmission_amplitude, beam_splitter(transmission_amplitude))
end

struct Rotation <: Interferometer
    angle::Float64
    U::Matrix
    Rotation(angle) = new(angle, rotation_matrix(angle))
end

struct PhaseShift <: Interferometer
    shifted_modes::Array
    param_::Array
    m::Int
    U::Matrix
    PhaseShift(shifted_modes, param_) = new(shifted_modes, param_, length(shifted_modes), phase_shift(shifted_modes, param_))
end
