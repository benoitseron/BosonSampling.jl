# Example usage: Generate random Gram matrix of given rank as input state

using Revise
using BosonSampling
using LinearAlgebra
using ArgCheck

# Parameters
n = 2  # number of photons
m = 3  # number of modes
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

# Let's construct an interferometer W
# It will be of size r*m and split each photon 1...n into groups
# Each group dof_basis = 1...r represents the components of photons 1...n on the basis vector indexed by dof_basis
# Each group has m modes. The top n < m modes are filled, one by one. Mode i correspond to the coefficient of photon i in the dof_basis

r = r_effective

######### to be checked!

###### TODO add a check to see where the photons initially are. the code below assumes that the are in modes 1...n

@argcheck is_input_in_first_modes(input)

W = zeros(ComplexF64, m, r*m) # m input modes with n photons, output separates in r groups of m modes 

### conventions reminder ###

#   3. Scattering matrix M:
#   index_input = fill_arrangement(input_state)
#   index_output = fill_arrangement(output_state)
#   M = U[index_input, index_output]  # n×n submatrix
#     - Extracts rows corresponding to input photons
#     - Extracts columns corresponding to output photons

for l in 1:r # partial distinguishability basis 
    for k in 1:n # photon index 
        W[k, (l-1)*m+k] = V[k, l]
    end
end

W

V

#### TODO but now W is not unitary... 