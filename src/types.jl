include("special_matrices.jl")

abstract type Interferometer end

is_unitary(Interferometer) = true

struct UserDefinedInterferometer <: Interferometer
    m::Int
    U::Matrix
    UserDefinedInterferometer(U::Matrix) = is_unitary(U) ? new(m,U) : error("input matrix is not unitary")
end

struct RandHaar <: Interferometer
    m::Int
    U::Matrix{ComplexF64}
    RandHaar(m) = new(m, rand_haar(m))
end

struct Fourier <: Interferometer
    m::Int
    U::Matrix{ComplexF64}
    Fourier(m) = new(m, fourier_matrix(m))
end

struct BeamSplitter <: Interferometer
    m::Int
    reflectivity::Float64
    U::Matrix
    BeamSplitter(m,reflectivity) = new(m, reflectivity, beam_splitter(m, reflectivity))
end

struct Rotation <: Interferometer
    m::Int
    angle::Float64
    U::Matrix
    Rotation(m,angle) = new(m, angle, rotation_matrix(m, angle))
end

struct PhaseShift <: Interferometer
    m::Int
    param::Float64
    U::Matrix
    PhaseShift(m, param) = new(m, param, phase_shift(m, param))
end

struct ModeOccupation
    m::Int
    n::Int
    state::Vector{Int}
    ModeOccupation(state) = all(state[:] .>= 0) ? new(length(state), sum(state), state) : error("negative photon counts")
end

at_most_one_photon_per_bin(r::ModeOccupation) = all(r.state[:] .<= 1)
first_modes(m::Int, n::Int) = ModeOccupation([i <= n ? 1 : 0 for i in 1:m])

abstract type InputType end

abstract type Bosonic <: InputType end

abstract type PartDist <: InputType end

abstract type ToyModel <: PartDist end

abstract type RandomModel <: PartDist end

abstract type Classical <: InputType end

abstract type Undef <: InputType end

struct GramMatrix{T<:InputType}

    n::Int
    S::Matrix

    function GramMatrix{T}(n::Int) where {T<:InputType}

        if T == Bosonic
            return new{T}(n, ones(ComplexF64, n, n))
        elseif T == Classical
            return new{T}(n, Matrix{ComplexF64}(I, n, n))
        elseif T == Undef
            return new{T}(n, Matrix{ComplexF64}(undef, n, n))
        elseif T == RandomModel
            return new{T}(n, rand_gram_matrix(n))
        else
            error("type ", T, " not implemented")
        end

    end

    function GramMatrix{T}(n::Int, S::Matrix) where {T<:InputType}

        if T <: PartDist && T != RandomModel
            return new{T}(n, S)
        else
            T in [Bosonic, Classical, Undef, RandomModel] ? error("S matrix should not be specified for type ", T) : error("Type ", T, " not implemented")
        end

    end
end

struct Input{T<:InputType}

    r::ModeOccupation
    G::GramMatrix

    function Input{T}(r::ModeOccupation) where {T<:InputType}

        if T in [Bosonic, Classical, Undef, RandomModel]
            return new{T}(r, GramMatrix{T}(r.n))
        else
            error("Type ", T, " not implemented")
        end

    end

    function Input{T}(r::ModeOccupation, G::GramMatrix) where {T<:InputType}

        if T <: PartDist && T != RandomModel
            return new{T}(r, G)
        else
            error("Type ", T, " not implemented")
        end

    end

end

abstract type OutputMeasurementType end

struct FockDetection <: OutputMeasurementType end

struct OutputMeasurement{T<:OutputMeasurementType}

    # as in most of boson sampling literature, this is for detectors blind to
    # internal degrees of freedom

    s::ModeOccupation # position of the detectors

    function OutputMeasurement{T}(s) where {T<:OutputMeasurementType}
        if T == FockDetection
            return at_most_one_photon_per_bin(s) ? new(s) : error("more than one detector per more")
        else
            return error(T, " not implemented")
        end
    end

end

mutable struct EventProbability
    probability::Union{Number,Nothing}
    precision::Union{Number,Nothing}
    # precision is set to eps() when doing non-randomised methods
    # although it is of course larger and this should be implemented
    # with permanent approximations, see for instance https://arxiv.org/abs/1904.06229
    failure_probability::Union{Number,Nothing}

    function EventProbability()
        new(nothing,nothing,nothing)
    end

    function EventProbability(probability)

        try
            probability = clean_proba(probability)
            new(probability,nothing,nothing)
        catch
            error("invalid probability")
        end
    end
end


struct Event{TIn<:InputType, TOut<:OutputMeasurementType}

    input_state::Input{TIn}
    output_measurement::OutputMeasurement{TOut}
    proba_params::EventProbability
    interferometer::Interferometer

    function Event{TIn, TOut}(input_state, output_measurement, interferometer::Interferometer, proba_params = nothing) where {TIn<:InputType, TOut<:OutputMeasurementType}

        if proba_params == nothing
            proba_params = EventProbability()
        end

        new{TIn,TOut}(input_state, output_measurement, proba_params, interferometer)
    end

    Event(i, o, interf, p=nothing) = Event{get_parametric_type(i)[1], get_parametric_type(o)[1]}(i, o, interf, p)

end

# i = Input{Bosonic}(first_modes(3,4))
# o = OutputMeasurement{FockDetection}(first_modes(3,4))
# ev = Event(i,o, EventProbability())
#
# @show ev.proba_params.probability
