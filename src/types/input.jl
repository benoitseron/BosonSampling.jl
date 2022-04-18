abstract type InputType end

struct Bosonic <: InputType
end

abstract type PartDist <: InputType end
struct OneParameterInterpolation <: PartDist
end
struct RandomGramMatrix <: PartDist
end
struct UserDefinedGramMatrix <: PartDist
end

struct Distinguishable <: InputType
end

struct Undef <: InputType
end

abstract type Gaussian <: InputType end

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
    distinguishability_param::Union{Float64,Nothing}
    generating_vectors::OrthonormalBasis

    function GramMatrix{T}(n::Int) where {T<:InputType}
        if T == Bosonic
            return new{T}(n, ones(ComplexF64,n,n), nothing, nothing, OrthonormalBasis())
        elseif T == Distinguishable
            return new{T}(n, Matrix{ComplexF64}(I,n,n), nothing, nothing, OrthonormalBasis())
        elseif T == RandomGramMatrix
            return new{T}(n, rand_gram_matrix(n), nothing, nothing, OrthonormalBasis())
        elseif T == Undef
            return new{T}(n, Matrix{ComplexF64}(undef,n,n), nothing, nothing, OrthonormalBasis())
        else
            error("type ", T, " not implemented")
        end
    end

    function GramMatrix{T}(n::Int, distinguishability_param::Float64) where {T<:InputType}
        if T == OneParameterInterpolation
            return new{T}(n, gram_matrix_one_param(n,distinguishability_param), nothing, distinguishability_param, OrthonormalBasis())
        else
            T in [Bosonic,Distinguishable,RandomGramMatrix,Undef] ? error("S matrix should not be specified for type ", T) : error("type ", T, " not implemented")
        end
    end

    function GramMatrix{T}(n::Int, S::Matrix) where {T<:InputType}
        if T == UserDefinedGramMatrix
            return new{T}(n, S, nothing, nothing, OrthonormalBasis())
        else
            T in [Bosonic,Distinguishable,RandomGramMatrix,Undef] ? error("S matrix should not be specified for type ", T) : error("type ", T, " not implemented")
        end
    end
    
end

struct Input{T<:InputType}

    r::ModeOccupation
    n::Int
    m::Int
    G::GramMatrix
    distinguishability_param::Union{Float64,Nothing}

    function Input{T}(r::ModeOccupation) where {T<:InputType}
        if T in [Bosonic, Distinguishable, Undef, RandomGramMatrix]
            return new{T}(r, r.n, r.m, GramMatrix{T}(r.n), nothing)
        else
            error("type ", T, " not implemented")
        end
    end

    function Input{T}(r::ModeOccupation, distinguishability_param::Float64) where {T<:InputType}
        if T == OneParameterInterpolation
            return new{T}(r, r.n, r.m, GramMatrix{T}(r.n,distinguishability_param), distinguishability_param)
        else
            error("type ", T, " not implemented")
        end
    end

    function Input{T}(r::ModeOccupation, S::Matrix) where {T<:InputType}
        if T == UserDefinedGramMatrix
            return new{T}(r, r.n, r.m, GramMatrix{T}(r.n,S), nothing)
        else
            error("type ", T, " not implemented")
        end
    end

end
