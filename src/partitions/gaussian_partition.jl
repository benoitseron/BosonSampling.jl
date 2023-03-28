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

input_state = GeneralGaussian(m = m)
interferometer = RandHaar(m)
part = equilibrated_partition(m, 1)

# part = Partition(Subset(first_modes(1, m)))
# part.subsets[1].subset

function compute_probabilities_partition_gaussian(physical_interferometer::Interferometer, part::Partition, input_state::GeneralGaussian)


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


        ### cutoff ###

        n_max = 10 #div(2m, part.n_subset) 
        @warn "arbitrary cutoff"


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

                Q = zeros(ComplexF64, 4 .* size(C))

                Q[1:m, 1:m] = I - C 
                Q[1:m, 2m+1:3m] = - C - virtual_interferometer_matrix
                Q[1:m, 3m+1:4m] = -1im .* virtual_interferometer_matrix

                Q[m+1:2m, m+1:2m] = I + C 
                Q[m+1:2m, 2m+1:3m] = -1im .* virtual_interferometer_matrix
                Q[m+1:2m, 3m+1:4m] = - C + virtual_interferometer_matrix

                Q[2m+1:3m, 1:m] = - C - virtual_interferometer_matrix
                Q[2m+1:3m, m+1:2m] = -1im .* virtual_interferometer_matrix
                Q[2m+1:3m, 2m+1:3m] = I - C

                Q[3m+1:4m, 1:m] = -1im .* virtual_interferometer_matrix
                Q[3m+1:4m, m+1:2m] = - C + virtual_interferometer_matrix
                Q[3m+1:4m, 3m+1:4m] = I + C


                coeffs = prod([2 / sqrt((1+ λ[j] * exp(2*r[j]))*(1+ λ[j] * exp(-2*r[j]))) for j in 1:m])

                probas_fourier[index_fourier_array] = coeffs * (det(Q))^(-0.5) * exp(0.5 * dot(Λ, Q^(-1) * Λ) - dot(Λ_plus, delta_x) - dot(Λ_minus, delta_y))
                
        end

        physical_indexes = copy(fourier_indexes)

        probas_physical(physical_index) = 1/(n_max+1)^(part.n_subset) * sum(probas_fourier[i] * exp(-2pi*1im/(n_max+1) * dot(physical_index, fourier_index)) for (i,fourier_index) in enumerate(fourier_indexes))

        pdf = [probas_physical(physical_index) for physical_index in physical_indexes]

        pdf = clean_pdf(pdf)

        mc = MultipleCounts(ModeOccupation.(physical_indexes), pdf)

        mc

end

mc = compute_probabilities_partition_gaussian(interferometer, part, input_state)
pdf = mc.proba

bar(real(pdf))

sort_samples_total_photon_number_in_partition!(mc)

bar(mc.proba)

### fft ###

# using FFTW

# probas_fourier

# probas_fourier_matrix = reshape(probas_fourier, ((n_max + 1), div(length(probas_fourier), (n_max + 1))))

# pdf_matrix = ifft(probas_fourier_matrix)

# pdf = reshape(pdf_matrix, (length(probas_fourier),))

# clean_pdf(pdf)

# bar(real.(pdf), alpha = 0.5)

# a = randn(20)

# ifft(fft(a))


### test case with one mode ###

# check if n is even

"""

        pdf_one_mode(n,r)

Returns the probability of having n photons in one mode, given the displacement r. Input of 1 mode squeezed sttate.
"""
pdf_one_mode(n,r::Real) = begin
        if n % 2 == 1
                return 0.0
        else

                μ = cosh(r)
                ν = sinh(r)

                return 1/μ * (ν/μ)^(n) * factorial(n) / (2^(div(n,2)) * factorial(div(n,2)))^2
        end
end

pdf_one_mode_array = [pdf_one_mode(i,r[1]) for i in 0:n_max]

bar!(real(pdf_one_mode_array), alpha = 0.5)

