"""
Corrected Clifford & Clifford Boson Sampling Algorithm Implementation

Based on the theoretical framework from arXiv:1706.01260 and correct boson sampling principles.
Key insights:
- O(n 2^n + poly(m,n)) complexity
- Sequential sampling approach
- Proper permanent calculations for conditional probabilities
"""

using LinearAlgebra
using StatsBase
using Random

"""
    clifford_algorithm_corrected(U::Matrix, input_modes::Vector{Int})

Corrected implementation of the Clifford & Clifford algorithm.
U: m×m unitary matrix
input_modes: vector of length n indicating which modes have input photons
"""
function clifford_algorithm_corrected(U::Matrix, input_modes::Vector{Int})
    m = size(U, 1)
    n = length(input_modes)
    
    # Extract relevant submatrix - columns for input modes
    A = U[:, input_modes]
    
    # Randomly permute input photons (this is important for the algorithm)
    perm = randperm(n)
    A_perm = A[:, perm]
    
    # Initialize output mode sequence
    output_modes = Int[]
    
    # Sample first output mode: probability ∝ |A[i,1]|²
    probs_1 = abs2.(A_perm[:, 1])
    mode_1 = wsample(1:m, Weights(probs_1))
    push!(output_modes, mode_1)
    
    # Sample remaining output modes sequentially
    for k in 2:n
        # For the k-th photon, compute conditional probabilities
        probs_k = zeros(Float64, m)
        
        # The probability for mode i is proportional to:
        # |sum_{j=1}^k A[i,j] * Minor_{j}|²
        # where Minor_j is the permanent of the (k-1)×(k-1) submatrix
        # obtained by removing column j and using rows from already sampled modes
        
        for i in 1:m
            amplitude = zero(ComplexF64)
            
            for j in 1:k
                # Compute the minor: permanent of submatrix without column j
                if k == 2
                    # For k=2, the "submatrix" is just a 1×1 matrix
                    minor = A_perm[output_modes[1], j == 1 ? 2 : 1]
                else
                    # For k>2, we need permanent of (k-1)×(k-1) submatrix
                    cols_to_keep = [1:(j-1); (j+1):k]
                    rows_to_use = output_modes[1:(k-1)]
                    
                    if length(cols_to_keep) > 0 && length(rows_to_use) > 0
                        submatrix = A_perm[rows_to_use, cols_to_keep]
                        minor = permanent_fast(submatrix)
                    else
                        minor = one(ComplexF64)
                    end
                end
                
                amplitude += A_perm[i, j] * minor
            end
            
            probs_k[i] = abs2(amplitude)
        end
        
        # Sample the k-th output mode
        if sum(probs_k) > 1e-12
            mode_k = wsample(1:m, Weights(probs_k))
        else
            # Fallback if all probabilities are zero (numerical issue)
            mode_k = rand(1:m)
        end
        push!(output_modes, mode_k)
    end
    
    # Convert mode sequence to occupancy vector
    occupancy = zeros(Int, m)
    for mode in output_modes
        occupancy[mode] += 1
    end
    
    return occupancy
end

"""
    permanent_fast(A::Matrix)

Fast permanent calculation using the appropriate method based on matrix size.
"""
function permanent_fast(A::Matrix)
    n = size(A, 1)
    
    if n == 0
        return one(eltype(A))
    elseif n == 1
        return A[1, 1]
    elseif n == 2
        return A[1,1]*A[2,2] + A[1,2]*A[2,1]
    else
        # Use the built-in permanent function for larger matrices
        return BosonSampling.permanent(A)
    end
end

"""
    clifford_sampler_corrected(input::Input, interf::Interferometer)

Corrected Clifford sampler for BosonSampling.jl types.
"""
function clifford_sampler_corrected(input::Input, interf::Interferometer)
    U = interf.U
    input_arrangement = BosonSampling.fill_arrangement(input)
    
    return clifford_algorithm_corrected(U, input_arrangement)
end

"""
    alternative_clifford_implementation(U::Matrix, input_modes::Vector{Int})

Alternative implementation based on the recursive structure of permanents.
This follows more closely the theoretical description from the literature.
"""
function alternative_clifford_implementation(U::Matrix, input_modes::Vector{Int})
    m = size(U, 1)
    n = length(input_modes)
    
    # Get the relevant submatrix
    A = U[:, input_modes]
    
    # Sample output modes one by one
    output_sequence = Int[]
    remaining_inputs = collect(1:n)
    
    for photon in 1:n
        # Compute probabilities for each output mode
        probs = zeros(Float64, m)
        
        for output_mode in 1:m
            # Probability is |permanent of augmented matrix|² / |permanent of current matrix|²
            
            if photon == 1
                # First photon: just |A[output_mode, input]|²
                prob_amplitude = A[output_mode, remaining_inputs[1]]
            else
                # For subsequent photons, we need to compute the conditional probability
                # This involves permanents of submatrices
                
                # Create the matrix for this configuration
                current_rows = [output_sequence; output_mode]
                current_cols = remaining_inputs[1:photon]
                
                if length(current_rows) == length(current_cols)
                    submatrix = A[current_rows, current_cols]
                    prob_amplitude = permanent_fast(submatrix)
                    
                    # Normalize by the permanent of the previous configuration
                    if photon > 1
                        prev_matrix = A[output_sequence, remaining_inputs[1:(photon-1)]]
                        prev_permanent = permanent_fast(prev_matrix)
                        if abs(prev_permanent) > 1e-12
                            prob_amplitude /= prev_permanent
                        end
                    end
                else
                    prob_amplitude = zero(ComplexF64)
                end
            end
            
            probs[output_mode] = abs2(prob_amplitude)
        end
        
        # Sample the output mode for this photon
        if sum(probs) > 1e-12
            sampled_mode = wsample(1:m, Weights(probs))
        else
            sampled_mode = rand(1:m)
        end
        
        push!(output_sequence, sampled_mode)
        
        # Remove the used input (this implements the random permutation implicitly)
        if photon < n
            remove_idx = rand(1:length(remaining_inputs))
            splice!(remaining_inputs, remove_idx)
        end
    end
    
    # Convert to occupancy vector
    occupancy = zeros(Int, m)
    for mode in output_sequence
        occupancy[mode] += 1
    end
    
    return occupancy
end

"""
    clifford_sampler_alternative(input::Input, interf::Interferometer)

Alternative implementation wrapper.
"""
function clifford_sampler_alternative(input::Input, interf::Interferometer)
    U = interf.U
    input_arrangement = BosonSampling.fill_arrangement(input)
    
    return alternative_clifford_implementation(U, input_arrangement)
end