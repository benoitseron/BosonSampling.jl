### Inputs ###

"""
Supertype to any concrete input type such as `Bosonic`
"""
abstract type InputType end

struct Bosonic <: InputType
end

"""
Type used to notify that the input is made of FockState partially distinguishable
photons.
"""
abstract type PartDist <: InputType
end

"""
One parameter model of partial distinguishability interpolating between indistinguishable
photons FockState and fully distinguishable.

see for the definition of the x-model :
[Sampling of partially distinguishable bosons and the relation to the
multidimensional permanent](https://arxiv.org/pdf/1410.7687.pdf)
"""
struct OneParameterInterpolation <: PartDist
end

"""
Model of partially distinguishable photons FockState described by a randomly generated Gram matrix.
"""
struct RandomGramMatrix <: PartDist
end

struct UserDefinedGramMatrix <: PartDist
end

struct Distinguishable <: InputType
end

struct Undef <: InputType
end

"""
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
Matrix of partial distinguishability. Will automatically generate the proper
matrix if given an input type defining a specific gram matrix such as `Bosonic`.
"""
struct GramMatrix{T<:InputType}

    n::Int
    S::Matrix
    rank::Union{Int,Nothing}
    distinguishability_param::Union{Real,Nothing}
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

    function GramMatrix{T}(n::Int, distinguishability_param::Real) where {T<:InputType}
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

"""

    Input{T<:InputType}
    Input{T}(r::ModeOccupation) where {T<:InputType}
    Input{T}(r::ModeOccupation, G::GramMatrix) where {T<:InputType}

Input state at the entrance of the interferometer.
"""
struct Input{T<:InputType}

    r::ModeOccupation
    n::Int
    m::Int
    G::GramMatrix
    distinguishability_param::Union{Real,Nothing}

    function Input{T}(r::ModeOccupation) where {T<:InputType}
        if T in [Bosonic, Distinguishable, Undef, RandomGramMatrix]
            return new{T}(r, r.n, r.m, GramMatrix{T}(r.n), nothing)
        else
            error("type ", T, " not implemented")
        end
    end

    function Input{T}(r::ModeOccupation, distinguishability_param::Real) where {T<:InputType}
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

at_most_one_photon_per_bin(input::Input) = at_most_one_photon_per_bin(input.r)
