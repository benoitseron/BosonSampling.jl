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
    using Parameters
    using UnPack
    using Optim


    using Permanents
    using Plots
    using Test
    using Combinatorics
    using Random
    using IterTools
    using Statistics
    using LinearAlgebra #removed so as to be able to use generic types such as BigFloats, can be put back if needed
    using PolynomialRoots
    using StatsBase
    #using JLD
    using CSV
    using DataFrames
    using Tables
    using Plots#; plotly() #plotly allows to make "dynamic" plots where you can point the mouse and see the values, they can also be saved like that but note that this makes them heavy (and html instead of png)
    using PrettyTables #to display large matrices
    using Roots
    using BenchmarkTools
    using Optim
    using ProgressMeter
    using ProgressBars
    using Parameters
    using ArgCheck
    using Distributions
    using Luxor
    using AutoHashEquals
    using LinearRegression
    using HypothesisTests
    using SimpleTraits
    using Parameters
    using UnPack
    using Dates
    using JLD
    using DelimitedFiles
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
    using Parameters
    using UnPack
    using Combinatorics
    using BosonSampling
    using Plots
    using ProgressMeter
    using Distributions
    using Random
    using Test
    using ArgCheck
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
    using Parameters
    using UnPack
    using Optim
    using DelimitedFiles
    using Dates
    using DelimitedFiles    

    using ITensors
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

### interferometer ###
m = 1
#physical_interferometer = Fourier(m)
#U = physical_interferometer.U
U = ones(ComplexF64, m,m)


### squeezing ###
r = 0.5 * ones(m)
λ = 1 * ones(m) # no thermal noise is 1
displacement = (0.0 + 0.0im) * ones(m)
delta_x = real.(displacement)
delta_y = imag.(displacement)

### useful matrices ###

C_array = [0.5 - 1/(1+ λ[j] * exp(2*r[j])) for j in 1:m]
C = diagm(C_array)

### Lambda matrices ###

Λ_plus = [2 * delta_x[j] /(1+ λ[j] * exp(2*r[j])) for j in 1:m]
Λ_minus = [2 * delta_y[j] /(1+ λ[j] * exp(-2*r[j])) for j in 1:m]

Λ = vcat(Λ_plus , Λ_minus, Λ_plus, - Λ_minus)

### partition ###

#part = equilibrated_partition(m, 2)

part = equilibrated_partition(m, 1)

### cutoff ###

n_max = 10
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

                @show fourier_element

                this_phase = exp(2*pi*1im/(n_max+1) * fourier_element)

                @show this_phase

                for j in 1:length(diag)

                        if part.subsets[i].subset[j] == 1

                                diag[j] = this_phase - 1 # previously multiply by phase
                                @show diag[j]

                        end

                end

        end

        @show diag 

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

# pdf = clean_pdf(pdf)

bar(real(pdf), alpha = 0.5)

pdf

### test case with one mode ###

# check if n is even


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