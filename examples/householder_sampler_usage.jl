"""
Example: Using the Householder-based Partial Distinguishability Sampler

This demonstrates how to use sample!() with UserDefinedGramMatrix inputs
for arbitrary partial distinguishability specified by Gram matrices.
"""

using BosonSampling
using LinearAlgebra

println("="^70)
println("Householder Partial Distinguishability Sampler - Usage Example")
println("="^70)

# Setup parameters
n = 3  # number of photons
m = 5  # number of modes
r = 2  # rank of Gram matrix (degree of distinguishability)

# Method 1: Random Gram matrix with specified rank
# ================================================
println("\n### Method 1: Random Gram matrix ###\n")

S_random = rand_gram_matrix_from_orthonormal_basis(n, r)
input_random = Input{UserDefinedGramMatrix}(first_modes(n, m), S_random)
interf = RandHaar(m)
ev_random = Event(input_random, FockSample(), interf)

sample!(ev_random)
println("Sampled output: ", ev_random.output_measurement.s)


# Method 2: Custom Gram matrix
# =============================
println("\n### Method 2: Custom Gram matrix ###\n")

# Create a custom Gram matrix (must be positive semi-definite, Hermitian)
# Example: Two photons highly indistinguishable, third one more distinguishable
S_custom = ComplexF64[
    1.0      0.95     0.3;
    0.95     1.0      0.35;
    0.3      0.35     1.0
]

println("Custom Gram matrix:")
display(S_custom)
println()

input_custom = Input{UserDefinedGramMatrix}(first_modes(n, m), S_custom)
ev_custom = Event(input_custom, FockSample(), interf)

sample!(ev_custom)
println("\nSampled output: ", ev_custom.output_measurement.s)


# Method 3: Comparing different distinguishability levels
# ========================================================
println("\n### Method 3: Comparing distinguishability levels ###\n")

# Fully indistinguishable (rank 1)
S_indist = ones(ComplexF64, n, n)
input_indist = Input{UserDefinedGramMatrix}(first_modes(n, m), S_indist)

# Partially distinguishable (rank 2)
S_partial = rand_gram_matrix_from_orthonormal_basis(n, 2)
input_partial = Input{UserDefinedGramMatrix}(first_modes(n, m), S_partial)

# Fully distinguishable (rank n)
S_dist = Matrix{ComplexF64}(I, n, n)
input_dist = Input{UserDefinedGramMatrix}(first_modes(n, m), S_dist)

# Sample from each
samples_indist = []
samples_partial = []
samples_dist = []

for i in 1:5
    ev = Event(input_indist, FockSample(), interf)
    sample!(ev)
    push!(samples_indist, ev.output_measurement.s.state)

    ev = Event(input_partial, FockSample(), interf)
    sample!(ev)
    push!(samples_partial, ev.output_measurement.s.state)

    ev = Event(input_dist, FockSample(), interf)
    sample!(ev)
    push!(samples_dist, ev.output_measurement.s.state)
end

println("Fully indistinguishable samples:")
for (i, s) in enumerate(samples_indist)
    println("  Sample $i: $s")
end

println("\nPartially distinguishable samples:")
for (i, s) in enumerate(samples_partial)
    println("  Sample $i: $s")
end

println("\nFully distinguishable samples:")
for (i, s) in enumerate(samples_dist)
    println("  Sample $i: $s")
end

println("\n" * "="^70)
println("Key Points:")
println("="^70)
println("  • Use Input{UserDefinedGramMatrix}(occupation, S) with your Gram matrix S")
println("  • S[i,j] = ⟨photon_i, photon_j⟩ represents photon overlaps")
println("  • S must be positive semi-definite and Hermitian")
println("  • rank(S) determines the degree of distinguishability")
println("  • Use sample!() as with any other input type")
println("="^70)
