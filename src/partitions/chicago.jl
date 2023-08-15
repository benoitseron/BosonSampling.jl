include("gaussian_partition.jl")

# using the formalism of Quantum-inspired classical algorithm for molecular vibronic spectra




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