include("gaussian_partition.jl")

# using the formalism of Quantum-inspired classical algorithm for molecular vibronic spectra




function compute_probabilities_partition_gaussian_chicago(physical_interferometer::Interferometer, part::Partition, input_state::GeneralGaussian, n_max = 11)


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

    # mc

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

m = 1
input_state = GeneralGaussian(m = m, r = 0.6 * ones(m)) 
interferometer = RandHaar(m)
part = equilibrated_partition(m, 1)
n_max = 50

mc = compute_probabilities_partition_gaussian_chicago(interferometer, part, input_state, n_max)

bar(real.(mc.proba))

sum(mc.proba)



@testset "gaussian partition - chicago algorithm" begin
    
    @testset "single mode" begin
            
            
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

                        return convert(Float64,1/μ * (ν/μ)^(n) * factorial(big(n)) / (2^(div(n,2)) * factorial(big(div(n,2))))^2)
                end
        end

        m = 1
        input_state = GeneralGaussian(m = m, r = 0.3 * ones(m))
        interferometer = Fourier(m)
        part = equilibrated_partition(m, 1)

        # part = Partition(Subset(first_modes(1, m)))
        # part.subsets[1].subset

        n_max = 51
        mc = compute_probabilities_partition_gaussian_chicago(interferometer, part, input_state, n_max)

        pdf = mc.proba

        bar(pdf, alpha = 0.5)

        pdf_one_mode_array = [pdf_one_mode(i,input_state.r[1]) for i in 0:n_max]

        # bar(pdf, alpha = 0.5)
        # bar!(pdf_one_mode_array, alpha = 0.5)

        if length(pdf) != length(pdf_one_mode_array)
            error("lengths are not equal")
        else

            #bar(pdf .- pdf_one_mode_array)

            @test pdf ≈ pdf_one_mode_array atol = 1e-10
        end

    end

end

