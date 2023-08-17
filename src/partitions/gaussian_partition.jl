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


"""

    convergence_checks_gaussian_partition(mc, probas_fourier)

Verifies that there is no error in the computation of a Gaussian binned counts distribution.
"""
function convergence_checks_gaussian_partition(mc, probas_fourier)

    sort_samples_total_photon_number_in_partition!(mc) # important for symmetry checks

    if any(real.(mc.proba) .< -ATOL) 

        # for proba in mc.proba
        #     if real.(proba) .< -ATOL
        #         println("negative probability: $proba")
        #     end
        # end 
        
        error("negative probabilities")
    
    end

    @argcheck isapprox(sum(mc.proba), 1., atol = ATOL) "not normalized"
    
    ### fourier symmetry checks ###

    first_half_range = 1: div(length(probas_fourier),2) + 1
    imag_fourier_first_half = (imag.(fftshift(probas_fourier)))[first_half_range]

    # problems seem to arise when the imaginary part is not asymmetric - this property means that the real space coefficients have a non zero imaginary part

    symmetric_part = 0.5 * (imag_fourier_first_half .+ reverse(imag_fourier_first_half))

    @argcheck isapprox(symmetric_part, zeros(length(symmetric_part)), atol = ATOL) "imaginary probabilies - the imaginary part of the Fourier transform is not asymmetric"

end


function compute_probabilities_partition_gaussian(physical_interferometer::Interferometer, part::Partition, input_state::GeneralGaussian, n_max = 11)

        if n_max % 2 == 0
                @warn "n_max must be odd for FFT purposes, converting"
                n_max = n_max +1
        end

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

        convergence_checks_gaussian_partition(mc, probas_fourier)

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

    find_cutoff(input_state
    include("gaussian_partition.jl")::GeneralGaussian; atol = ATOL)

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

### plotting the average number of photons distribution ###

"""
    function partition_expectation_values_gaussian(input_state::GeneralGaussian, n_max::Int, part::Partition)

Computes the asymptotic law of the expectation values of a Haar averaged partition for Gaussian input states .
"""
function partition_expectation_values_gaussian(input_state::GeneralGaussian, n_max::Int, part::Partition)


    #generating the list of indexes for counts 
    if n_max % 2 == 0
        @warn "n_max must be odd for FFT purposes, converting"
        n_max = n_max +1
    end

    counts = all_mode_configurations(n_max,part.n_subset, only_photon_number_conserving = false)

    partition_occupancy_array = [PartitionOccupancy(count, part) for count in counts]

    # need to only compute the expectation values for even total photon number because of input gaussian states

    is_total_photon_number_even(partition_occupancy::PartitionOccupancy) = sum(partition_occupancy.counts.state) % 2 == 0

    partition_expectation_values_array_dist = [is_total_photon_number_even(occ) ? probability_n_photons(occ.n, input_state) * partition_expectation_values(occ)[1] : 0 for occ in partition_occupancy_array]

    partition_expectation_values_array_bosonic = [is_total_photon_number_even(occ) ? probability_n_photons(occ.n, input_state) * partition_expectation_values(occ)[2] : 0 for occ in partition_occupancy_array]

    mc_dist =  MultipleCounts(ModeOccupation.(counts), partition_expectation_values_array_dist)
    mc_bosonic =  MultipleCounts(ModeOccupation.(counts), partition_expectation_values_array_bosonic)

    sort_samples_total_photon_number_in_partition!(mc_dist)
    sort_samples_total_photon_number_in_partition!(mc_bosonic)

    mc_bosonic, mc_dist

end



function compute_probabilities_partition_gaussian_chicago(physical_interferometer::Interferometer, part::Partition, input_state::GeneralGaussian, n_max = 11; give_debug_info = false)


        if !all(input_state.displacement .== 0.0)
            error("displacement not implemented")
        end
    
    
        if n_max % 2 == 0
            @warn "n_max must be odd for FFT purposes, converting"
            n_max = n_max +1
        end
    
        physical_interferometer = interferometer
        ### function ###
    
        @warn "arbitrary cutoff"
    
        # unpack the parameters of the input state
    
        @unpack m, r, λ, delta_x, delta_y = input_state
    
    
        renormalization_factor = 4^m 
        @warn "normalization factor should be one - to debug!"
    
    
    
        γ = @. [-1 + exp(2*r[j]) for j in 1:m]
    
        U = physical_interferometer.U
    
        ### useful matrices ###
    
        Γ = diagm(γ)
        N = prod(@. sqrt(1+ γ) / (π * γ))
    
        ### virtual interferometer matrix ###
    
        fourier_indexes = all_mode_configurations(n_max,part.n_subset, only_photon_number_conserving = false)
        probas_fourier = Array{ComplexF64}(undef, length(fourier_indexes))
        top_right_matrix = similar(U)
        bottom_left_matrix = similar(U)
    
        Q_non_psd = []
    
        for (index_fourier_array, fourier_index) in enumerate(fourier_indexes)
    
            # for each fourier index, we recompute the virtual interferometer
            top_right_matrix  = conj(U)
            bottom_left_matrix = U
    
            diag = [0.0 + 0im for i in 1:m] # previously 1.0 + 0im
            
            for (i,fourier_element) in enumerate(fourier_index)
    
                    # @show fourier_element
    
                    this_phase = exp(2*pi*1im/(n_max+1) * fourier_element)
    
                    # @show this_phase
    
                    for j in 1:length(diag)
                            # @show i,j
                            if part.subsets[i].subset[j] == 1
    
                                    diag[j] = this_phase 
    
                            end
    
                    end
    
            end
    
            # @show diag 
    
            top_right_matrix *= Diagonal(diag)
            bottom_left_matrix *= Diagonal(diag)
    
            top_right_matrix *= - conj(U') # more practical than the transpose function 
            bottom_left_matrix *= - U' ################ check minus sign
    
            ### matrix Q ###
    
            Q = zeros(ComplexF64, 2m, 2m)
    
            Q[1:m, 1:m] = 2 * Γ^(-1) + diagm([1.0 + 0im for i in 1:m])
            Q[1:m, m+1:2m] = top_right_matrix
            Q[m+1:2m, 1:m] = bottom_left_matrix
            Q[m+1:2m, m+1:2m] = 2 * Γ^(-1) + diagm([1.0 + 0im for i in 1:m])
    
            if !is_positive_semidefinite(real.(Q))
                @warn "real part of Q is not positive semidefinite for fourier index $fourier_index"
                push!(Q_non_psd, Q)
            end
    
            probas_fourier[index_fourier_array] = N * (2\pi)^(m) * det(Q)^(-0.5)
            
        end
    
        if Q_non_psd != []
            @warn "the real part of Q is not positive semidefinite for some fourier indexes"
    
            for Q in Q_non_psd
                @show Q
                @show eigvals(Q)
            end
    
            return nothing
        end
    
        physical_indexes = copy(fourier_indexes)
    
        probas_fourier_matrix = reshape(probas_fourier, ((n_max + 1), div(length(probas_fourier), (n_max + 1))))
    
        shifted_probas_fourier_matrix = fftshift(probas_fourier_matrix)
        pdf_matrix = fft(shifted_probas_fourier_matrix) / (n_max + 1)^part.n_subset
    
        pdf = renormalization_factor * reshape(pdf_matrix, (length(probas_fourier),))
    
        try
            pdf = clean_pdf(pdf)
        catch
            @warn "invalid pdf, skipping cleaning"
            # println("press any key to continue")
            # readline()
            pdf = real.(pdf)
        end
    
        mc = MultipleCounts(ModeOccupation.(physical_indexes), pdf)

        convergence_checks_gaussian_partition(mc, probas_fourier)
    
        if give_debug_info
            return mc, probas_fourier, pdf_matrix, shifted_probas_fourier_matrix, pdf
        else
            return mc
        end
    
    end