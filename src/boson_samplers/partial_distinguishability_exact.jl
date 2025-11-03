# Example usage: Generate random Gram matrix of given rank as input state

using Revise
using BosonSampling
using LinearAlgebra

# Parameters
n = 3  # number of photons
m = 5  # number of modes
r = 2  # rank of Gram matrix (r < n for partial distinguishability)

# Generate random Gram matrix of rank r
# This represents the overlap matrix between n photons in r-dimensional space
S = rand_gram_matrix_from_orthonormal_basis(n, r)

# Create input state with partial distinguishability
input = Input{UserDefinedGramMatrix}(first_modes(n, m), S)

# Create interferometer
interf = RandHaar(m)

# Extract internal degrees of freedom from Gram matrix
# S[i,j] = ⟨photon_i, photon_j⟩ = Σₖ V[i,k] * conj(V[j,k])
# where V[i,k] is the k-th coefficient of photon i in the internal basis
V = gram_to_coefficients(S)

r_effective = size(V, 2)

# Verify: reconstruct Gram matrix from V
S_reconstructed = reconstruct_gram_matrix(V)

# Consistency Check:
# Original S and Reconstructed S = V * V' should be equal
@assert S ≈ S_reconstructed "Gram matrix reconstruction failed!"

println("\n✓ Successfully extracted $r_effective internal degrees of freedom from Gram matrix")
println("✓ Gram matrix reconstruction verified")

