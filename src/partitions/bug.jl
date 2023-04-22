include("gaussian_partition.jl")

### does not fail ###
begin

    m = 20
    input_state = GeneralGaussian(m = m, r = 0.40 * ones(m))
    interferometer = RandHaar(m)
    part = equilibrated_partition(m, 2)

    # part = Partition(Subset(first_modes(1, m)))
    # part.subsets[1].subset

    n_max = 50
    mc = compute_probabilities_partition_gaussian(interferometer, part, input_state, n_max)

    # bar(real(pdf))

    sort_samples_total_photon_number_in_partition!(mc)

    @show mc

    bar(mc.proba)

end


### fails ###
begin

    m = 20
    input_state = GeneralGaussian(m = m, r = 0.42 * ones(m))
    interferometer = RandHaar(m)
    part = equilibrated_partition(m, 2)

    # part = Partition(Subset(first_modes(1, m)))
    # part.subsets[1].subset

    n_max = 25
    mc = compute_probabilities_partition_gaussian(interferometer, part, input_state, n_max)

    # bar(real(pdf))

    sort_samples_total_photon_number_in_partition!(mc)

    @show mc

    bar(mc.proba)

end

bar(mc.proba)

sum(mc.proba[500:end])

bar(mc.proba[300:end])

### single mode ###

# also failures for single mode

# interesting behaviour: failure at 0.42, n_max = 50
# but not at 0.42, n_max = 30

# failure for n_max = 33, 32
# not failure n_max = 31


begin

    m = 20
    input_state = GeneralGaussian(m = m, r = 0.42 * ones(m))
    interferometer = RandHaar(m)
    part = equilibrated_partition(m, 1)

    # part = Partition(Subset(first_modes(1, m)))
    # part.subsets[1].subset

    n_max = 30
    mc = compute_probabilities_partition_gaussian(interferometer, part, input_state, n_max)

    # bar(real(pdf))

    sort_samples_total_photon_number_in_partition!(mc)

    @show mc

    bar(mc.proba)

end

### other trials ###


begin

    m = 1
    input_state = GeneralGaussian(m = m, r = 1.5 * ones(m))
    interferometer = RandHaar(m)
    part = equilibrated_partition(m, 1)

    # part = Partition(Subset(first_modes(1, m)))
    # part.subsets[1].subset

    n_max = 50
    mc = compute_probabilities_partition_gaussian(interferometer, part, input_state, n_max)

    # bar(real(pdf))

    sort_samples_total_photon_number_in_partition!(mc)

    @show mc

    bar(mc.proba)

end


###### zero padding ######

m = 20
input_state = GeneralGaussian(m = m, r = 0.42 * ones(m))
interferometer = RandHaar(m)
part = equilibrated_partition(m, 1)


n_cutoff = 30

# double n_max for padding purposes - we will not compute the unnecessary coefficients however

@warn "need to feed n_cutoff now!"

n_max = 2 * n_cutoff

if n_max % 2 == 0
    @warn "n_max must be odd for FFT purposes, converting"
    n_max = n_max + 1
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

begin
    fourier_indexes = all_mode_configurations(n_max,part.n_subset, only_photon_number_conserving = false)
    probas_fourier = Array{ComplexF64}(undef, length(fourier_indexes))
    virtual_interferometer_matrix = similar(U)


    for (index_fourier_array, fourier_index) in enumerate(fourier_indexes)

        # do not make the computations if more than the cutoff

        # if all(fourier_index .<= n_cutoff)

            println("computing for fourier index $(fourier_index)")

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

        # else

        #     # the zero padding for the fourier transform
        #     # probas_fourier[index_fourier_array] = 0.0 + 0im

        # end
        
    end
end
#probas_fourier[n_cutoff + 1 : end] .= 0.0 + 0im

probas_fourier

physical_indexes = copy(fourier_indexes)

probas_fourier_matrix = reshape(probas_fourier, ((n_max + 1), div(length(probas_fourier), (n_max + 1))))

shifted_probas_fourier_matrix = fftshift(probas_fourier_matrix)
pdf_matrix = fft(shifted_probas_fourier_matrix) / (n_max + 1)^part.n_subset

pdf = reshape(pdf_matrix, (length(probas_fourier),))

bar(real.(pdf), alpha = 0.5)
bar!(imag.(pdf) , alpha = 0.5)
bar!(abs.(pdf) , alpha = 0.5)

# probas_fourier_matrix = reshape(probas_fourier, ((n_max + 1), div(length(probas_fourier), (n_max + 1))))

# probas_fourier_matrix = zeros(eltype(probas_fourier_matrix_padded), size(probas_fourier_matrix_padded))

# padded_length = Int(ceil(3(n_max + 1) / 2))

# probas_fourier_padded = zeros(eltype(probas_fourier), (padded_length, div(length(probas_fourier), (n_max + 1))))

# probas_fourier_padded[1:(n_max + 1), :] = probas_fourier_matrix

shifted_probas_fourier_matrix = similar(probas_fourier_matrix_padded)

# shifted_probas_fourier_matrix[1:div(size(probas_fourier_matrix_padded)[1], 2), :] = fftshift(probas_fourier_matrix_padded[1:div(size(probas_fourier_matrix_padded)[1], 2), :])

shifted_probas_fourier_matrix = fftshift(probas_fourier_matrix_padded)

# shifted_probas_fourier_matrix = probas_fourier_matrix_padded

# pdf_matrix = fft(shifted_probas_fourier_matrix) / (n_max + 1)^part.n_subset

pdf_matrix = fft(shifted_probas_fourier_matrix) / (n_max + 1)^part.n_subset

pdf_matrix 

bar(real.(pdf_matrix))

pdf = reshape(pdf_matrix, (length(probas_fourier),))

bar(real.(pdf))

try
    pdf = clean_pdf(pdf)
catch
    @warn "invalid pdf, skipping cleaning"
    # println("press any key to continue")
    # readline()
    pdf = real.(pdf)
end

mc = MultipleCounts(ModeOccupation.(physical_indexes), pdf)

mcfunction compute_probabilities_partition_gaussian(physical_interferometer::Interferometer, part::Partition, input_state::GeneralGaussian, n_max = 11, padding_factor = 2)

if n_max % 2 == 0
    @warn "n_max must be odd for FFT purposes, converting"
    n_max = n_max + 1
end

...

fourier_indexes = all_mode_configurations(n_max, part.n_subset, only_photon_number_conserving = false)
probas_fourier = Array{ComplexF64}(undef, length(fourier_indexes))
virtual_interferometer_matrix = similar(U)

...

probas_fourier_matrix = reshape(probas_fourier, ((n_max + 1), div(length(probas_fourier), (n_max + 1))))

# Create a new matrix with the desired padded size
padded_probas_fourier_matrix = zeros(ComplexF64, (n_max + 1) * padding_factor, div(length(probas_fourier), (n_max + 1)))

# Copy the elements from the original matrix to the center of the new matrix
half_padding = div((n_max + 1) * (padding_factor - 1), 2)
padded_probas_fourier_matrix[half_padding + 1 : half_padding + n_max + 1, :] = probas_fourier_matrix

shifted_padded_probas_fourier_matrix = fftshift(padded_probas_fourier_matrix)
pdf_matrix = fft(shifted_padded_probas_fourier_matrix) / ((n_max + 1) * padding_factor)^part.n_subset

...

mc = MultipleCounts(ModeOccupation.(physical_indexes), pdf)

mc

end


### gpt4 ###

padding_factor = 2

# function compute_probabilities_partition_gaussian(physical_interferometer::Interferometer, part::Partition, input_state::GeneralGaussian, n_max = 11, padding_factor = 2)

    if n_max % 2 == 0
        @warn "n_max must be odd for FFT purposes, converting"
        n_max = n_max + 1
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

    fourier_indexes = all_mode_configurations(n_max, part.n_subset, only_photon_number_conserving = false)
    probas_fourier = Array{ComplexF64}(undef, length(fourier_indexes))
    virtual_interferometer_matrix = similar(U)

    for (index_fourier_array, fourier_index) in enumerate(fourier_indexes)

        # for each fourier index, we recompute the virtual interferometer
        virtual_interferometer_matrix  = conj(U)

        diag = [0.0 + 0im for i in 1:m] # previously 1.0 + 0im

        for (i, fourier_element) in enumerate(fourier_index)

            this_phase = exp(2*pi*1im/(n_max+1) * fourier_element)

            for j in 1:length(diag)
                if part.subsets[i].subset[j] == 1
                    diag[j] = this_phase - 1
                end
            end
        end

        virtual_interferometer_matrix *= Diagonal(diag)
        virtual_interferometer_matrix *= conj(U')

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

    probas_fourier_matrix = reshape(probas_fourier, ((n_max + 1) * padding_factor, div(length(probas_fourier), (n_max + 1))))

    shifted_probas_fourier_matrix = fftshift(probas_fourier_matrix)
    pdf_matrix = ifft(shifted_probas_fourier_matrix) / ((n_max + 1) * padding_factor)^part.n_subset

    pdf = reshape(pdf_matrix, (length(probas_fourier),))

    try
        pdf = clean_pdf(pdf)
    catch
        @warn "invalid pdf, skipping cleaning"
        pdf = real.(pdf)
    end

    mc = MultipleCounts(ModeOccupation.(physical_indexes), pdf)

    mc

# end





