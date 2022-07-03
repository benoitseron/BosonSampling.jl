### Inputs ###

"""
Supertype to any concrete input type such as [`Bosonic`](@ref), [`PartDist`](@ref), [`Distinguishable`](@ref)
and [`Undef`](@ref).
"""
abstract type InputType end

"""
Type used to notify that the input is made of FockState indistiguishable photons.
"""
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
and fully distinguishable photons FockState.

!!! note "Reference"
    [Sampling of partially distinguishable bosons and the relation to the
    multidimensional permanent](https://arxiv.org/pdf/1410.7687.pdf)
"""
struct OneParameterInterpolation <: PartDist
end

"""
Model of partially distinguishable photons FockState described by a randomly generated [`GramMatrix`](@ref).
"""
struct RandomGramMatrix <: PartDist
end

"""
Model of partially distinguishable photons FockState described by a provided [`GramMatrix`](@ref).
"""
struct UserDefinedGramMatrix <: PartDist
end

"""
Model of distinguishable photons FockState.
"""
struct Distinguishable <: InputType
end

"""
Model of photons FockState with undefined [`GramMatrix`](@ref).
"""
struct Undef <: InputType
end

"""
Type used to notify that the input is made of Gaussian states.
"""
abstract type Gaussian end

"""
Type used to notify that the input is made of the vacuum state.
"""
struct VacuumState <: Gaussian

    displacement::Vector{Complex}
    covariance_matrix::Matrix{Complex}

    function VacuumState()
        new([0;0], [1/2 0; 0 1/2])
    end

end

"""
Type used to notify that the input is made of coherent state.
"""
struct CoherentState <: Gaussian

    trunc::Int
    displacement_parameter::Complex
    displacement::Vector{Complex}
    covariance_matrix::Matrix{Complex}
    spectrum::Vector{Real}

    function CoherentState(trunc::Int, displacement::Complex)
        new(trunc,
            displacement_parameter,
            sqrt(2) * [displacement_parameter; conj(displacement_parameter)],
            [1/2 0; 0 1/2],
            [exp(-1/2*abs(displacement_parameter)^2) * displacement_parameter^j / sqrt(factorial(j)) for j in 0:trunc])
    end

end

"""
Type used to notify that the input is made of thermal state.
"""
struct ThermalState <: Gaussian

    trunc::Int
    mean_photon_number::Real
    displacement::Vector{Complex}
    covariance_matrix::Matrix{Complex}
    spectrum::Vector{Real}

    function ThermalState(trunc::Int, mean_photon_number::Real)
        new(trunc,
            0,
            mean_photon_number,
            zeros(2),
            [mean_photon_number+1/2 0; 0 mean_photon_number+1/2],
            [mean_photon_number^j / (1+mean_photon_number)^(j+1) for j in 0:trunc])
    end

end

"""
Type used to notify that the input is made of single mode squeezed state.
"""
struct SingleModeSqueezedVacuum <: Gaussian

    trunc::Int
    squeezing_parameter::Vector{Real}
    displacement::Vector{Complex}
    covariance_matrix::Matrix{Complex}
    spectrum::Vector{Real}

    function SingleModeSqueezedVacuum(trunc::Int, squeezing_parameter::Vector{Real})

        r = squeezing_parameter[1]
        θ = squeezing_parameter[2]
        spectrum = zeros(trunc+1)'
        for j in 1:length(spectrum)
            iseven(j) ? spectrum[j] = sqrt(sech(r)) * sqrt(factorial(2j))/factorial(j) * (-1/2*exp(θ*im)*tanh(r))^j : nothing
        end

        new(trunc,
            squeezing_parameter,
            zeros(2),
            1/2 * cosh(2r) * Matrix(I,2,2) - 1/2 * sinh(2r) * [cos(θ) sin(θ); sin(θ) -cos(θ)],
            spectrum)
    end

end

# struct TwoModeSqueezedVacuum <: Gaussian
#
#     trunc::Int
#     squeezing_parameter::Vector{Real}
#     displacement::Vector{Complex}
#     covariance_matrix::Matrix{Complex}
#     spectrum::Vector{Real}
#
#     function TwoModeSqueezedVacuum(trunc, squeezing_parameter)

"""
    OrthonormalBasis(vector_matrix::Union{Matrix, Nothing})

Basis of vectors ``v_1,...,v_n`` stored as columns in a ``n``-by-``r`` matrix
possibly empty.

    Fields:
        vectors_matrix::Union{Matrix,Nothing}
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
    GramMatrix{T}(n::Int) where {T<:InputType}
    GramMatrix{T}(n::Int, distinguishability_param::Real) where {T<:InputType}
    GramMatrix{T}(n::Int, S::Matrix) where {T<:InputType}

Matrix of partial distinguishability. Will automatically generate the proper
matrix related to the provided [`InputType`](@ref).

    Fields:
        - n::Int: photons number
        - S::Matrix: Gram matrix
        - rank::Union{Int, Nothing}
        - distinguishability_param::Union{Real, Nothing}
        - generating_vectors::OrthonormalBasis
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

    Fields:
        - r::ModeOccupation
        - n::Int
        - m::Int: modes numbers
        - G::GramMatrix
        - distinguishability_param::Union{Real, Nothing}
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
