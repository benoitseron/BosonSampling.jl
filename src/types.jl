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


struct ModeOccupation
    n::Int
    m::Int
    state::Vector{Int}
    ModeOccupation(state) = all(state[:] .>= 0) ? new(sum(state), length(state), state) : error("negative photon counts")
end

at_most_one_photon_per_bin(state) = all(state[:] .<= 1)
at_most_one_photon_per_bin(r::ModeOccupation) = at_most_one_photon_per_bin(r.state)
isa_subset(subset_modes::ModeOccupation) = at_most_one_photon_per_bin(subset_modes)
isa_subset(state) = at_most_one_photon_per_bin(state)
# this last function is defined to avoid the abuse of language
first_modes(n::Int,m::Int) = n<=m ? ModeOccupation([i <= n ? 1 : 0 for i in 1:m]) : error("n>m")

abstract type InputType end

struct Bosonic <:InputType
end
struct PartDist <:InputType
end
struct Distinguishable <:InputType
end
struct Undef <:InputType
end

mutable struct OrthonormalBasis
    """basis of vectors v_1,...,v_n stored as columns in a n*r matrix
    possibly empty"""
    vectors_matrix::Union{Matrix,Nothing}
    function OrthonormalBasis(vectors_matrix = nothing)
        if vectors_matrix == nothing
            new(nothing)
        else
            is_orthonormal(vectors_matrix, atol = ATOL) ? new(vectors_matrix) : error("invalid orthonormal basis")
        end
    end
end

struct GramMatrix{T<:InputType}

    n::Int
    S::Matrix
    rank::Union{Int,Nothing}
    generating_vectors::OrthonormalBasis

    function GramMatrix{T}(n::Int) where {T<:InputType}
        if T == Bosonic
            return new{T}(n,ones(ComplexF64,n,n), nothing, OrthonormalBasis())
        elseif T == Distinguishable
            return new{T}(n,Matrix{ComplexF64}(I,n,n), nothing, OrthonormalBasis())
        elseif T == Undef
            return new{T}(n,Matrix{ComplexF64}(undef,n,n), nothing, OrthonormalBasis())
        else
            error("type ", T, " not implemented")
        end
    end
    function GramMatrix{T}(n::Int,S::Matrix) where {T<:InputType}
        if T in [Bosonic, Distinguishable, Undef]
            error("S matrix should not be specified for [Bosonic, Distinguishable, Undef] types")
            ########### to be made better
        elseif T == PartDist
            return new{T}(n,S, nothing, OrthonormalBasis())
        else
            error("type ", T, " not implemented")
        end
    end
end

struct Input{T<:InputType}
    r::ModeOccupation
    G::GramMatrix
    n::Int
    m::Int
    function Input{T}(r::ModeOccupation) where {T<:InputType}
        if T in [Bosonic, Distinguishable, Undef]
            return new{T}(r,GramMatrix{T}(r.n),r.n,r.m)
        else
            error("type ", T, " not implemented")
        end
    end
    function Input{T}(r::ModeOccupation, G::GramMatrix) where {T<:InputType}

        if T == PartDist
            return new{T}(r,G,r.n,r.m)
        else
            error("type ", T, " not implemented")
        end
    end
end

at_most_one_photon_per_bin(inp::Input) = at_most_one_photon_per_bin(inp.r)

abstract type OutputMeasurementType end

struct FockDetection <: OutputMeasurementType
end

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
    precision::Union{Number,Nothing} # see remarks in conventions
    failure_probability::Union{Number,Nothing}


    function EventProbability(probability = nothing)

        if probability == nothing
            new(nothing, nothing, nothing)
        else
            try
                probability = clean_proba(probability)
                new(probability,nothing,nothing)
            catch
                error("invalid probability")
            end
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

    Event(i,o,interf,p = nothing) = Event{get_parametric_type(i)[1], get_parametric_type(o)[1]}(i,o,interf,p)

end
