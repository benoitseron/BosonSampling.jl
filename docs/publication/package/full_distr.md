```julia
# Set the model of partial distinguishability
T = OneParameterInterpolation

# Define an input of 3 photons placed
# among 5 modes
x = 0.74
in = Input{T}(first_modes(3,5), x)

# Interferometer
l = 0.63
U = RandHaar(in.m)

# Compute the full output statistics
p_theo, p_approx, p_samp
= noisy_distribution(input=in,loss=l,interf=U)
```