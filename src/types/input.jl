### Inputs ###

abstract type InputType end

abstract type Bosonic <: InputType end

abstract type GaussianState <: Bosonic end

struct ThermalState <: GaussianState

    mode_number::Int
    mean_photon_number::Float64
    mean_vector::Vector
    covariance_matrix::Matrix

    function ThermalState(mode_number::Int, mean_photon_number::Float64)
        return new(
            mode_number,
            mean_photon_number,
            zeros(mode_number),
            (mean_photon_number + 1 / 2) *
            Matrix{ComplexF64}(I, 2*mode_number, 2*mode_number),
        )
    end

end

abstract type PartDist <: InputType end

struct ToyModel <: PartDist
    distinguishability::Float64
    ToyModel(distinguishability::Float64) = new(distinguishability)
end

abstract type RandomModel <: PartDist end

abstract type Undef <: InputType end

abstract type Distinguishable <: InputType end

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

mutable struct GramMatrix{T<:InputType}

    n::Int
    S::Matrix
    rank::Union{Int,Nothing}
    generating_vectors::OrthonormalBasis

    function GramMatrix{T}(n::Int) where {T<:InputType}
        if T <: Bosonic
            return new{T}(
                n,
                ones(ComplexF64, n, n),
                nothing,
                OrthonormalBasis(),
            )
        elseif T == RandomModel
            return new{T}(
                n,
                rand_gram_matrix(n),
                nothing,
                OrthonormalBasis(),
            )
        elseif T == Distinguishable
            return new{T}(
                n,
                Matrix{ComplexF64}(I, n, n),
                nothing,
                OrthonormalBasis(),
            )
        elseif T == Undef
            return new{T}(
                n,
                Matrix{ComplexF64}(undef, n, n),
                nothing,
                OrthonormalBasis(),
            )
        else
            error("type ", T, " not implemented")
        end
    end

    function GramMatrix{T}(n::Int, distinguishability::Float64) where {T<:InputType}
        if T == ToyModel
            return new{T}(
                n,
                gram_matrix_toy_model(n, distinguishability),
                nothing,
                OrthonormalBasis(),
            )
        else
            T in [Bosonic, Distinguishable, RandomModel, Undef] ?
            error("S matrix should not be specified for type ", T) :
            error("Type ", T, " not implemented")
        end
    end

end

struct Input{T<:InputType}

    r::ModeOccupation
    G::Matrix
    n::Int
    m::Int
    distinguishability::Any

    function Input{T}(r::ModeOccupation) where {T<:InputType}
        if T == Bosonic
            return new{T}(r, GramMatrix{T}(r.n).S, r.n, r.m, 1.0)
        elseif T == Distinguishable
            return new{T}(r, GramMatrix{T}(r.n).S, r.n, r.m, 0.0)
        elseif T in [RandomModel, Undef]
            return new{T}(r, GramMatrix{T}(r.n).S, r.n, r.m, nothing)
        else
            error("type ", T, " not implemented")
        end
    end

    function Input{T}(r::ModeOccupation, distinguishability::Float64) where {T<:InputType}
        if T == ToyModel
            return new{T}(
                r,
                GramMatrix{T}(r.n, distinguishability).S,
                r.n,
                r.m,
                distinguishability,
            )
        else
            error("type ", T, " not implemented")
        end
    end

end

at_most_one_photon_per_bin(inp::Input) = at_most_one_photon_per_bin(inp.r)
