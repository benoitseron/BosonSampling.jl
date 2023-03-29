### fft ###

using FFTW

### input parameters ###

m = 20
input_state = GeneralGaussian(m = m, r = 0.4 * ones(m))
interferometer = RandHaar(m)
part = equilibrated_partition(m, 2)

# part = Partition(Subset(first_modes(1, m)))
# part.subsets[1].subset

n_max = 201

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

bar(real.(pdf), alpha = 0.5)

@warn "skipping cleaning of probabilities and health checks"
try
        pdf = clean_pdf(pdf)
catch
        @warn "invalid pdf, skipping cleaning"
        println("press any key to continue")
        readline()
        pdf = real.(pdf)
end

mc = MultipleCounts(ModeOccupation.(physical_indexes), pdf)

mc


### FFT buisness ###

probas_fourier

probas_fourier_matrix = reshape(probas_fourier, ((n_max + 1), div(length(probas_fourier), (n_max + 1))))

shifted_probas_fourier_matrix = fftshift(probas_fourier_matrix)
pdf_matrix = fft(shifted_probas_fourier_matrix) / (n_max + 1)^part.n_subset

pdf_fft = reshape(pdf_matrix, (length(probas_fourier),))

# clean_pdf(pdf)

bar!(real.(pdf_fft), alpha = 0.5)

@test pdf ≈ pdf_fft

