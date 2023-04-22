begin
    using Revise

    using BosonSampling
    using Plots
    using ProgressMeter
    using Distributions
    using Random
    using Test
    using ArgCheck
    using StatsBase
    using ColorSchemes
    using Interpolations
    using Dierckx
    using LinearAlgebra
    using PrettyTables
    using LaTeXStrings
    using JLD
    using AutoHashEquals
    using LinearRegression
    
    using DataStructures
    using FFTW
    # include package using @with_kw
    using Parameters
    using Roots

end


# r - squeezing
# x - generatinc funtion
# η - vector of fourier coefficients
# m - modes
# Q - 
# Λ -
# δ -
# λ - thermal parameters
# V = U_round

# the lambdas are created page 24
# Q page 23
# C just above


"""

        mutable struct GeneralGaussian

A generalized input of a multimode gaussian state. 

        Fields:
                - m::Int = 1  
                - r::Vector{<:Real} = 0.5 * ones(m) # squeezing
                - λ::Vector{<:Real} = 1 * ones(m) # thermal noise, no thermal noise is 1
                - displacement::Vector{<:Complex} = (0.0 + 0.0im) * ones(m) 
                - delta_x::Vector{<:Real} = real.(displacement) # displacement real part
                - delta_y::Vector{<:Real} = imag.(displacement) # displacement imaginary part

"""
@with_kw mutable struct GeneralGaussian

        m::Int = 1  
        r::Vector{<:Real} = 0.5 * ones(m) # squeezing
        λ::Vector{<:Real} = 1 * ones(m) # thermal noise, no thermal noise is 1
        displacement::Vector{<:Complex} = (0.0 + 0.0im) * ones(m) 
        delta_x::Vector{<:Real} = real.(displacement) # displacement real part
        delta_y::Vector{<:Real} = imag.(displacement) # displacement imaginary part

end


function compute_probabilities_partition_gaussian(physical_interferometer::Interferometer, part::Partition, input_state::GeneralGaussian, n_max = 11)

        if n_max % 2 == 0
                @warn "n_max must be odd for FFT purposes, converting"
                n_max = n_max +1
        end

        physical_interferometer = interferometer
        ### function ###

        @warn "arbitrary cutoff"

        # unpack the parameters of the input state

        @unpack m, r, λ, delta_x, delta_y = input_state

        U = physical_interferometer.U

        ### useful matrices ###

        C_array = [0.5 - 1/(1+ λ[j] * exp(2*r[j])) for j in 1:m]
        C = diagm(C_array)

        ### Lambda matrices ###

        Λ_plus = [2 * delta_x[j] /(1+ λ[j] * exp(2*r[j])) for j in 1:m]
        Λ_minus = [2 * delta_y[j] /(1+ λ[j] * exp(-2*r[j])) for j in 1:m]

        Λ = vcat(Λ_plus , Λ_minus, Λ_plus, - Λ_minus)

        ### virtual interferometer matrix ###

        fourier_indexes = all_mode_configurations(n_max,part.n_subset, only_photon_number_conserving = false)
        probas_fourier = Array{ComplexF64}(undef, length(fourier_indexes))
        virtual_interferometer_matrix = similar(U)


        for (index_fourier_array, fourier_index) in enumerate(fourier_indexes)

                # for each fourier index, we recompute the virtual interferometer
                virtual_interferometer_matrix  = conj(U)

                diag = [0.0 + 0im for i in 1:m] # previously 1.0 + 0im
                
                for (i,fourier_element) in enumerate(fourier_index)

                        # @show fourier_element

                        this_phase = exp(2*pi*1im/(n_max+1) * fourier_element)

                        # @show this_phase

                        for j in 1:length(diag)
                                # @show i,j
                                if part.subsets[i].subset[j] == 1

                                        diag[j] = this_phase - 1 # previously multiply by phase
                                        # @show diag[j]

                                end

                        end

                end

                # @show diag 

                virtual_interferometer_matrix *= Diagonal(diag)
                virtual_interferometer_matrix *= conj(U') # more practical than the transpose function 

                ### matrix Q ###

                virtual_interferometer_matrix_transpose = conj(virtual_interferometer_matrix')

                Q = zeros(ComplexF64, 4 .* size(C))

                Q[1:m, 1:m] = I - C 
                Q[1:m, 2m+1:3m] = - C - virtual_interferometer_matrix
                Q[1:m, 3m+1:4m] = -1im .* virtual_interferometer_matrix

                Q[m+1:2m, m+1:2m] = I + C 
                Q[m+1:2m, 2m+1:3m] = -1im .* virtual_interferometer_matrix
                Q[m+1:2m, 3m+1:4m] = - C + virtual_interferometer_matrix

                Q[2m+1:3m, 1:m] = - C - virtual_interferometer_matrix_transpose
                Q[2m+1:3m, m+1:2m] = -1im .* virtual_interferometer_matrix_transpose
                Q[2m+1:3m, 2m+1:3m] = I - C

                Q[3m+1:4m, 1:m] = -1im .* virtual_interferometer_matrix_transpose
                Q[3m+1:4m, m+1:2m] = - C + virtual_interferometer_matrix_transpose
                Q[3m+1:4m, 3m+1:4m] = I + C


                coeffs = prod([2 / sqrt((1+ λ[j] * exp(2*r[j]))*(1+ λ[j] * exp(-2*r[j]))) for j in 1:m])

                probas_fourier[index_fourier_array] = coeffs * (det(Q))^(-0.5) * exp(0.5 * dot(Λ, Q^(-1) * Λ) - dot(Λ_plus, delta_x) - dot(Λ_minus, delta_y))
                
        end

        physical_indexes = copy(fourier_indexes)

        probas_fourier_matrix = reshape(probas_fourier, ((n_max + 1), div(length(probas_fourier), (n_max + 1))))

        shifted_probas_fourier_matrix = fftshift(probas_fourier_matrix)
        pdf_matrix = fft(shifted_probas_fourier_matrix) / (n_max + 1)^part.n_subset

        pdf = reshape(pdf_matrix, (length(probas_fourier),))

        try
                pdf = clean_pdf(pdf)
        catch
                @warn "invalid pdf, skipping cleaning"
                # println("press any key to continue")
                # readline()
                pdf = real.(pdf)
        end

        mc = MultipleCounts(ModeOccupation.(physical_indexes), pdf)

        mc

end

function probability_n_photons(n::Int, input_state::GeneralGaussian)

        m = input_state.m
    
        @argcheck m % 2 == 0 "only implemented for even m, still holds but need to use gamma functions instead of factorials"
        @argcheck all(input_state.r .== input_state.r[1]) "r must be the same for all modes"
        r = input_state.r[1]
        
        p(k, m) = binomial(div(m, 2) + k - 1, k) * sech(r)^m * tanh(r)^(2k)
    
        if n % 2 == 0
            return p(div(n,2), m)
        else
            return 0.
        end
    
end


"""

    find_cutoff(input_state::GeneralGaussian; atol = ATOL)

Given a value `atol`, finds the number of photons that is less likely to be detected than `atol`. Follows S3 in 2110.0156. 
"""
function find_cutoff(input_state::GeneralGaussian; atol = ATOL)

    n_in = input_state.m
    # check that all the elements of r are the same
    @argcheck all(input_state.r .== input_state.r[1]) "r must be the same for all modes"
    r = input_state.r[1]

    # @show n_in,r
    cutoff(α, n_in::Int, r::Real) = α * n_in * sinh(r)^2

    bound(α) = exp(-0.5 *α * n_in * (1-1/α)^2)

    f = x -> bound(x) - atol

    α = find_zero(f, 1.) 

    2*Int(ceil(cutoff(α, n_in, r))) # multiplied by two because this is for photon pairs

end