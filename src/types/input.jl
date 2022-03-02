### Inputs ###

abstract type InputType end

struct Bosonic <:InputType
end
struct PartDist <:InputType
end
abstract type ToyModel <: PartDist
end
abstract type RandomModel <: PartDist
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
        elseif T == RandomModel
            return new{T}(n, rand_gram_matrix(n), nothing, OrthonormalBasis())
        elseif T == Undef
            return new{T}(n,Matrix{ComplexF64}(undef,n,n), nothing, OrthonormalBasis())
        else
            error("type ", T, " not implemented")
        end
    end
    function GramMatrix{T}(n::Int,S::Matrix) where {T<:InputType}
        if T <: PartDist && T != RandomModel
            return new{T}(n, S, nothing, OrthonormalBasis())
        else
            T in [Bosonic, Distinguishable, Undef, RandomModel] ? error("S matrix should not be specified for type ", T) : error("Type ", T, " not implemented")
        end
    end
end

struct Input{T<:InputType}
    r::ModeOccupation
    G::GramMatrix
    n::Int
    m::Int
    function Input{T}(r::ModeOccupation) where {T<:InputType}
        if T in [Bosonic, Distinguishable, Undef, RandomModel]
            return new{T}(r,GramMatrix{T}(r.n),r.n,r.m)
        else
            error("type ", T, " not implemented")
        end
    end
    function Input{T}(r::ModeOccupation, G::GramMatrix) where {T<:InputType}

        if T <: PartDist && T != RandomModel
            return new{T}(r,G,r.n,r.m)
        else
            error("type ", T, " not implemented")
        end
    end
end

at_most_one_photon_per_bin(inp::Input) = at_most_one_photon_per_bin(inp.r)
