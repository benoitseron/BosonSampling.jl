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
Supertype to any concrete input type of Gaussian state.
"""
abstract type Gaussian end

"""
    VacuumState()

Type used to notify that the input is made of the vacuum state.

    Fields:
        displacement::Vector{Complex}
        covariance_matrix::Matrix{Complex}
"""
struct VacuumState <: Gaussian

    displacement::Vector{Float64}
    covariance_matrix::Matrix{Float64}

    function VacuumState()
        new(zeros(2), [1.0 0.0; 0.0 1.0])
    end

end

"""
    CoherentState(displacement_parameter::Complex)

Type used to notify that the input is made of a coherent state.

    Fields:
        - displacement_parameter::Complex
        - displacement::Vector{Complex}
        - covariance_matrix::Matrix{Complex}
"""
struct CoherentState <: Gaussian

    displacement_parameter::Complex
    displacement::Vector{Float64}
    covariance_matrix::Matrix{Float64}

    function CoherentState(displacement_parameter::Complex)
        new(displacement_parameter,
            sqrt(2) * [real(displacement_parameter); imag(displacement_parameter)],
            [1/2 0; 0 1/2])
    end

end

"""
    ThermalState(mean_photon_number::Real)

Type used to notify that the input is made of a thermal state.

    Fields:
        mean_photon_number::Real
        displacement::Vector{Complex}
        covariance_matrix::Matrix{Complex}
"""
struct ThermalState <: Gaussian

    mean_photon_number::Real
    displacement::Vector{Float64}
    covariance_matrix::Matrix{Float64}

    function ThermalState(mean_photon_number::Real)
        new(mean_photon_number,
            zeros(2),
            [mean_photon_number+1/2 0; 0 mean_photon_number+1/2])
    end

end

"""
    SingleModeSqueezedVacuum(squeezing_parameter::Real)

Type used to notify that the input is made of single mode squeezed state.

    Fields:
        - squeezing_parameter::Real
        - displacement::Vector{Complex}
        - covariance_matrix::Matrix{Complex}
"""
struct SingleModeSqueezedVacuum <: Gaussian

    squeezing_parameter::Real
    displacement::Vector{Float64}
    covariance_matrix::Matrix{Float64}

    function SingleModeSqueezedVacuum(squeezing_parameter::Real)

        new(squeezing_parameter,
            zeros(2),
            [exp(-2*squeezing_parameter) 0;
            0 exp(2*squeezing_parameter)])
    end

end

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

    function Input{T}(r::ModeOccupation, n::Int, m::Int, G::GramMatrix, distinguishability_param::Union{Real,Nothing}) where {T<:InputType}
        new{T}(r,n,m,G, distinguishability_param)
    end

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

"""
    GaussianInput{T<:Gaussian}
    GaussianInput{T}(r::ModeOccupation, squeezing_parameters::Vector, source_transmission::Union{Vector, Nothing})

Input state made off Gaussian states at the entrance of the interferometer.
Models of partial distinguishability for Gausian states being of different nature than those for Fock states,
we distinct input of [`Bosonic`](@ref) Fock states with indistinguishable Gaussian states.

!!! note
    The notion of [`ModeOccupation`](@ref) here is also different. Even though
    the object is the same, the mode occuputation here as to be understood as a "boolean"
    where the value `0` tells us that the mode is fed with the vacuum while `1`
    states that the mode contains a Gaussian state of type `T`.
"""
struct GaussianInput{T<:Gaussian}

    r::ModeOccupation
    n::Int
    m::Int
    displacement::Vector{Float64}
    covariance_matrix::Matrix{Float64}

    displacement_parameters::Union{Vector, Nothing}
    mean_photon_numbers::Union{Vector, Nothing}
    squeezing_parameters::Union{Vector, Nothing}

    function GaussianInput{T}() where {T<:Gaussian}
        if T == VacuumState
            return Input{Bosonic}(first_modes(0,r.m))
        else
            error("type ", T, " not implemented")
        end
    end

    function GaussianInput{CoherentState}(r::ModeOccupation, displacement_parameters::Vector)

        r.state[1] == 1 ? sq = CoherentState(displacement_parameters[1]) : sq = VacuumState()
        cov = sq.covariance_matrix
        R = sq.displacement

        for i in 2:r.m
            if r.state[i] == 1
                sq = CoherentState(displacement_parameters[i])
                R = [R;sq.displacement]
                    sigma = sq.covariance_matrix
            else
                vacc = VacuumState()
                R = [R;vacc.displacement]
                sigma = vacc.covariance_matrix
            end
            cov = direct_sum(cov, sigma)
        end

        return new{CoherentState}(r, r.n, r.m, R, cov, displacement_parameters, nothing, nothing, nothing)

    end


    # function GaussianInput{T}(r::ModeOccupation, mean_photon_numbers::Vector) where {T<:Gaussian}
    #
    #     if T == ThermalState
    #         r.state[1] == 1 ? sq = ThermalState(mean_photon_numbers[1]) : sq = VacuumState()
    #         cov = sq.covariance_matrix
    #         R = sq.displacement
    #
    #         for i in 2:r.m
    #             if r.state[i] == 1
    #                 sq = ThermalState[mean_photon_numbers[i]]
    #                 R = [R;sq.displacement]
    #                 sigma = sq.covariance_matrix
    #             else
    #                 vacc = VacuumState()
    #                 R = [R;vacc.displacement]
    #                 sigma = vacc.covariance_matrix
    #             end
    #             cov = direct_sum(cov, sigma)
    #         end
    #         return new{T}(r, r.n, r.m, R, cov, nothing, displacement_parameters, nothing, nothing)
    #     else
    #         error("type ", T, " not implemented")
    #     end
    #
    # end


    function GaussianInput{SingleModeSqueezedVacuum}(r::ModeOccupation, squeezing_parameters::Vector)

        r.state[1] == 1 ? sq = SingleModeSqueezedVacuum(squeezing_parameters[1]) : sq = VacuumState()
        cov = sq.covariance_matrix
        R = sq.displacement

        for i in 2:r.m
            if r.state[i] == 1
                sq = SingleModeSqueezedVacuum(squeezing_parameters[i])
                R = [R;sq.displacement]
                sigma = sq.covariance_matrix
            else
                vacc = VacuumState()
                R = [R;vacc.displacement]
                sigma = vacc.covariance_matrix
            end
            cov = direct_sum(cov, sigma)
        end
        return new{SingleModeSqueezedVacuum}(r, r.n, r.m, R, cov, nothing, nothing, squeezing_parameters)
    end

end

"""
    get_spectrum(state::Gaussian, k::Int)

Get the spectrum in the Fock basis of a `state` up to `k` photons.
"""
function get_spectrum(state::Gaussian, k::Int)

    if typeof(state) == VacuumState
        return 1

    elseif typeof(state) == CoherentState
        α = state.displacement_parameter
        return [exp(-0.5*abs(α)^2) * (α^n)/sqrt(factorial(n)) for n in 0:k]

    elseif typeof(state) == ThermalState
        μ = state.mean_photon_number
        return [μ^j / (1+μ)^(j+1) for j in 0:k]

    elseif typeof(state) == SingleModeSqueezedVacuum
        r = state.squeezing_parameter
        res = Vector{Float64}(undef, k)
        for i in 0:k
            if iseven(i)
                res[i] = sqrt(sech(r)) * sqrt(factorial(2n))/factorial(n) * (-0.5*tanh(r))^n
            else
                res[i] = 0
            end
        end
        return res
    end

end
