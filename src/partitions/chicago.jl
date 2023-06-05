include("gaussian_partition.jl")

# using the formalism of Quantum-inspired classical algorithm for molecular vibronic spectra

m = 10
input_state = GeneralGaussian(m = m, r = 0.4 * ones(m))
interferometer = RandHaar(m)
part = equilibrated_partition(m, 1)
n_max = 25


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
    bottom_left_matrix *= - U'

    ### matrix Q ###

    Q = zeros(ComplexF64, 2m, 2m)

    Q[1:m, 1:m] = 2 * Γ^(-1) + diagm([1.0 + 0im for i in 1:m])
    Q[1:m, m+1:2m] = top_right_matrix
    Q[m+1:2m, 1:m] = bottom_left_matrix
    Q[m+1:2m, m+1:2m] = 2 * Γ^(-1) + diagm([1.0 + 0im for i in 1:m])

    if !is_positive_semidefinite(Q)
        @warn "Q is not positive semidefinite for fourier index $fourier_index"
        push!(Q_non_psd, Q)
    end

    probas_fourier[index_fourier_array] = N * (2\pi)^(m) * det(Q)^(-0.5)
    
end

Q = Q_non_psd[1]
eigvals(Q)

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

bar(real.(mc.proba))

sum(mc.proba)