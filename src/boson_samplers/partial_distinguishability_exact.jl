# Example usage: Generate random Gram matrix of given rank as input state

using Revise
using BosonSampling
using LinearAlgebra
using ArgCheck
using SparseArrays

# Parameters
n = 4  # number of photons
m = 5  # number of modes
r = 3  # rank of Gram matrix (r < n for partial distinguishability)

# Generate random Gram matrix of rank r
# This represents the overlap matrix between n photons in r-dimensional space
S = rand_gram_matrix_from_orthonormal_basis(n, r)

# Create input state with partial distinguishability
input = Input{UserDefinedGramMatrix}(first_modes(n, m), S)

# Create interferometer
interf = RandHaar(m)

interf.U

# Extract internal degrees of freedom from Gram matrix
# S[i,j] = ⟨photon_i, photon_j⟩ = Σₖ V[i,k] * conj(V[j,k])
# where V[i,k] is the k-th coefficient of photon i in the internal basis
C = gram_to_coefficients(S)

r_effective = size(C, 2)

# Verify: reconstruct Gram matrix from C
S_reconstructed = reconstruct_gram_matrix(C)

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
@argcheck m >= n

# Householder transformation for photon i: V_i = I - 2*C_i*C_i† (r×r matrix)
# where C_i is the i-th row of C (representing photon i's internal structure)
splitting_matrix(i) = Matrix{ComplexF64}(I, r, r) - 2 * C[i,:] * C[i,:]'

# Enlarge to m×m for compatibility with the m×m physical interferometer U
# Only the first r×r block is non-trivial (active internal DOF)
# The remaining (m-r) modes are just identity (inactive)
function splitting_matrix_enlarged(i)
    result = Matrix{ComplexF64}(I, m, m)
    result[1:r, 1:r] = splitting_matrix(i)
    return result
end

# Create block diagonal with n blocks (one per photon)
# Each block is the m×m Householder transformation for that photon
# Total size: n*m × n*m
function get_all_splitting_matrices()
    result = Matrix{ComplexF64}(I, n*m, n*m)

    for photon in 1:n
        range_ = 1+(photon-1)*m : m+(photon-1)*m
        result[range_, range_] = splitting_matrix_enlarged(photon)
    end
    result
end

splitting_interferometer = get_all_splitting_matrices()

# Mode Shuffle Permutation
# ========================
# Purpose: Reorder modes so that modes at the same position within each block are grouped together
#
# Current layout (after splitting matrices): n blocks, each with m modes (one block per photon)
#   [Block₁: modes 1...m | Block₂: modes 1...m | ... | Blockₙ: modes 1...m]
#   Total: n*m modes
#
# Desired layout: Group by mode position across all blocks
#   [Mode₁ of all blocks | Mode₂ of all blocks | ... | Modeₘ of all blocks]
#   = [B₁M₁, B₂M₁, ..., BₙM₁, B₁M₂, B₂M₂, ..., BₙM₂, ..., B₁Mₘ, B₂Mₘ, ..., BₙMₘ]
#
# This allows each group to be processed by a separate copy of interferometer U
#
# Indexing formula:
#   For a mode at linear position k ∈ [1, n*m]:
#   - Block index (photon):  block = floor((k-1)/m) + 1  ∈ [1, n]
#   - Mode in block:         mode = ((k-1) mod m) + 1    ∈ [1, m]
#   - Old position:          k_old = (block-1)*m + mode  (row-major: block varies slowest)
#   - New position:          k_new = (mode-1)*n + block  (column-major: block varies fastest)
#
# Example with n=2, m=3:
#   Old: [B₁M₁, B₁M₂, B₁M₃, B₂M₁, B₂M₂, B₂M₃] = [1, 2, 3, 4, 5, 6]
#   New: [B₁M₁, B₂M₁, B₁M₂, B₂M₂, B₁M₃, B₂M₃] = [1, 4, 2, 5, 3, 6]
#
# This is equivalent to reshape-transpose-reshape:
#   1. View n*m linear index as (n, m) matrix in row-major order
#   2. Transpose to (m, n)
#   3. Flatten in row-major order

function mode_shuffle_permutation()
    P = zeros(ComplexF64, n*m, n*m)
    for block in 1:n
        for mode in 1:m
            k_old = (block-1)*m + mode
            k_new = (mode-1)*n + block
            P[k_new, k_old] = 1.0
        end
    end
    return P
end

shuffle_permutation = mode_shuffle_permutation()


function block_diagonal_interferometers()
    result = Matrix{ComplexF64}(I, m*n, m*n)

    for pd_group in 1:n
        range_ = 1+(pd_group-1)*m: m+(pd_group-1)*m
        result[range_, range_] = interf.U
    end
    result
end

block_diagonal_interferometer = block_diagonal_interferometers()

full_interferometer = shuffle_permutation' * block_diagonal_interferometer* shuffle_permutation * splitting_interferometer
# note: we unshuffled the modes - the physical mode 1 corresponds to the first m modes summed over, etc for the rest 

occupation_simulated_photons = zeros(Int, m*n)

for photon in 1:n 
    occupation_simulated_photons[(photon-1)*m + 1] = 1  
end



input_simulated = Input{Bosonic}(ModeOccupation(occupation_simulated_photons))
interf_simulated = UserDefinedInterferometer(full_interferometer)

# Sample using Clifford's algorithm for indistinguishable bosons
# In the expanded space, photons are treated as bosonic
output_sample = FockSample()
ev_simulated = Event(input_simulated, output_sample, interf_simulated)
sample!(ev_simulated)

# Get the sampled output in the expanded space (n*m modes)
sampled_expanded = ev_simulated.output_measurement.s

# Convert ModeOccupation to vector for indexing
sampled_expanded_vec = sampled_expanded.state

# Binning: Collapse back to m physical modes
# Each physical mode i corresponds to modes {i, m+i, 2m+i, ..., (n-1)m+i} in expanded space
# Sum the photon counts from all n blocks for each mode position
sampled_physical = zeros(Int, m)

for photon_block in 1:n
    for mode in 1:m
        expanded_mode = (photon_block-1)*m + mode
        sampled_physical[mode] += sampled_expanded_vec[expanded_mode]
    end
end

sampled_physical

println("\nSampled output (expanded space): ", sampled_expanded)
println("Sampled output (physical space): ", sampled_physical)
println("Total photons (should be $n): ", sum(sampled_physical))

