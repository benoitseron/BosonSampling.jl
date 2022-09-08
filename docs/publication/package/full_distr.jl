# Set the model of partial distinguishability
T = OneParameterInterpolation

# Define an input of 3 photons placed
# among 5 modes
x = 0.74
i = Input{T}(first_modes(3,5), x)

# Interferometer
l = 0.63
U = RandHaar(i.m)

# Compute the full output statistics
p_exact, p_truncated, p_sampled =
noisy_distribution(input=i,loss=l,interf=U)
