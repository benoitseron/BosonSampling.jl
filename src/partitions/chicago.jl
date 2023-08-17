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

m = 20
input_state = GeneralGaussian(m = m, r = 0.6 * ones(m)) 
interferometer = RandHaar(m)
part = equilibrated_partition(m, 2)
n_max = 10

mc, probas_fourier = compute_probabilities_partition_gaussian_chicago(interferometer, part, input_state, n_max, give_debug_info = true)
