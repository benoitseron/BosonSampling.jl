### Inputs ###

"""
    InputType

Supertype to any concrete input type such as `Bosonic`
"""
abstract type InputType end

"""

    Bosonic

Type use to notify that the input is made of FockState indistinguishable photons.
"""
struct Bosonic <:InputType
end

"""

    PartDist

Type use to notify that the input is made of FockState partially distinguishable
photons.
"""
struct PartDist <:InputType
end
abstract type ToyModel <: InputType
end
abstract type RandomModel <: InputType
end
struct Distinguishable <:InputType
end
struct Undef <:InputType
end

"""

    OrthonormalBasis

Basis of vectors v_1,...,v_n stored as columns in a n*r matrix
possibly empty.
"""
mutable struct OrthonormalBasis

    vectors_matrix::Union{Matrix,Nothing}
    function OrthonormalBasis(vectors_matrix = nothing)
        if vectors_matrix == nothing
            new(nothing)
        else
            is_orthonormal(vectors_matrix, atol = ATOL) ? new(vectors_matrix) : error("invalid orthonormal basis")
        end
    end
end

"""

    GramMatrix{T<:InputType}

Matrix of partial distinguishability. Will automatically generate the proper
matrix if given an input type defining a specific gram matrix such as `Bosonic`.
"""
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

"""

    Input{T<:InputType}
    Input{T}(r::ModeOccupation) where {T<:InputType}
    Input{T}(r::ModeOccupation, G::GramMatrix) where {T<:InputType}

Input state at the entrance of the interferometer.
"""
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

        if T in [PartDist, ToyModel]
            return new{T}(r,G,r.n,r.m)
        else
            error("type ", T, " not implemented")
        end
    end
end

at_most_one_photon_per_bin(inp::Input) = at_most_one_photon_per_bin(inp.r)
