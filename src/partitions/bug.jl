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

    n_max = 15
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


### looking at inside of function ###

plt_failure = plot()


m = 20
input_state = GeneralGaussian(m = m, r = 0.42 * ones(m))
interferometer = RandHaar(m)
part = equilibrated_partition(m, 1)

n_max = 16


# begin

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

store_eigvals = []
store_Q = nothing

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

        # @warn "removing transpose for debugging purposes"   
        # virtual_interferometer_matrix_transpose = virtual_interferometer_matrix

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

        @show is_positive_semidefinite(Q)
        
        store_eigvals = eigvals(Q)
        store_Q = Q
        
    end

# end

store_Q

is_positive_semidefinite(store_Q)

fourier_indexes[3]

store_eigvals
eigenvalues = store_eigvals
real.(eigenvalues)

abs.(imag.(eigenvalues))


plt_failure = plot()

bar!(plt_failure, real.(fftshift(probas_fourier)), label = "real, n = $n_max", alpha = 0.5)
bar!(plt_failure, imag.(fftshift(probas_fourier)), label = "imag, n = $n_max" , alpha = 0.5)


plt_failure = plot()

bar!(plt_failure, abs.(fftshift(probas_fourier)), label = "n = $n_max", alpha = 0.5)
bar!(plt_failure, angle.(fftshift(probas_fourier)), label = "n = $n_max" , alpha = 0.5)

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

bar(mc.proba)




plt_failure