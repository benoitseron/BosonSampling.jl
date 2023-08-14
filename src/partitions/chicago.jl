include("gaussian_partition.jl")

# using the formalism of Quantum-inspired classical algorithm for molecular vibronic spectra




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

    if give_debug_info
        return mc, probas_fourier, pdf_matrix, shifted_probas_fourier_matrix, pdf
    else
        return mc
    end

end

#### tests #### 

# buggy case:
# m = 10
# input_state = GeneralGaussian(m = m, r = 0.65 * ones(m)) 
# interferometer = RandHaar(m)
# part = equilibrated_partition(m, 2)
# n_max = 25

# another buggy case

# m = 10
# input_state = GeneralGaussian(m = m, r = 0.7 * ones(m)) 
# interferometer = RandHaar(m)
# part = equilibrated_partition(m, 1)
# n_max = 100

m = 10
input_state = GeneralGaussian(m = m, r = 0.6 * ones(m)) 
interferometer = RandHaar(m)
part = equilibrated_partition(m, 2)
n_max = 20

mc, probas_fourier = compute_probabilities_partition_gaussian_chicago(interferometer, part, input_state, n_max, give_debug_info = true)

sort_samples_total_photon_number_in_partition!(mc)

bar(real.(mc.proba))

sum(mc.proba)

probas_fourier

begin
    
    plt_failure = plot()

    bar!(plt_failure, real.(fftshift(probas_fourier)), label = "real, n = $n_max", alpha = 0.5)
    bar!(plt_failure, imag.(fftshift(probas_fourier)), label = "imag, n = $n_max" , alpha = 0.5)
end


# given the properties of the FFT, one value is ommited in the array by construction
# this is why I take only the first half of the array
########### need to check carefully what we expect here 
########### also maybe the Fourier transform is twice too long?

first_half_range = 1: div(length(probas_fourier),2) + 1
imag_fourier_first_half = (imag.(fftshift(probas_fourier)))[first_half_range]

bar(imag_fourier_first_half)

# problems seem to arise when the imaginary part is not asymmetric - this property means that the real space coefficients have a non zero imaginary part

symmetric_part = 0.5 * (imag_fourier_first_half .+ reverse(imag_fourier_first_half))

@argcheck isapprox(symmetric_part, zeros(length(symmetric_part)), atol = ATOL) "the imaginary part of the Fourier transform is not asymmetric"



bar(symmetric_part)