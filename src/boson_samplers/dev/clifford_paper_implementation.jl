"""
# Clifford & Clifford Boson Sampling Algorithm - Paper Implementation

This is Algorithm B from the Clifford & Clifford paper (arXiv:1706.01260).
Notebook-style implementation for easy checking and understanding.

## Algorithm B: Boson sampler single sample z from q(z) in O(n2^n + poly(m,n)) time

**Require:** m and n positive integers; A first n columns of m × m Haar random unitary matrix

1. r ← ∅                                          ⊳ EMPTY ARRAY
2. A ← PERMUTE(A)                                  ⊳ RANDOMLY PERMUTE COLUMNS OF A  
3. wi ← |ai,1|², i ∈ [m]                          ⊳ MAKE INDEXED WEIGHT ARRAY w
4. x ← SAMPLE(w)                                   ⊳ SAMPLE INDEX x FROM w
5. r ← (r, x)                                      ⊳ APPEND x TO r
6. FOR k ← 2 TO n DO
7.     B_k^x ← A^[k]                              ⊳ Submatrix of A
8.     COMPUTE {Per B_k,ℓ^x, ℓ ∈ [k]}            ⊳ AS LEMMA 2
9.     wi ← |∑_{ℓ=1}^k ai,ℓ Per B_k,ℓ^x|², i ∈ [m]  ⊳ USING LAPLACE EXPANSION
10.    x ← SAMPLE(w)
11.    r ← (r, x)
12. END FOR
13. z ← INCSORT(r)                                ⊳ SORT r IN NON-DECREASING ORDER
14. RETURN z
"""

using LinearAlgebra
using StatsBase
using Random

#%% Cell 1: Helper Functions

"""
    permute_columns(A::Matrix, n::Int)

Randomly permute the first n columns of matrix A.
This corresponds to line 2: A ← PERMUTE(A)
"""
function permute_columns(A::Matrix, n::Int)
    m = size(A, 1)
    perm = randperm(n)
    A_permuted = zeros(eltype(A), m, n)
    for i in 1:n
        A_permuted[:, i] = A[:, perm[i]]
    end
    return A_permuted
end

"""
    compute_permanents_lemma2(B_k_x::Matrix, k::Int)

Compute {Per B_{k,ℓ}^x, ℓ ∈ [k]} as mentioned in line 8.
This computes permanents of submatrices with column ℓ removed.
"""
function compute_permanents_lemma2(B_k_x::Matrix, k::Int)
    permanents = Vector{ComplexF64}(undef, k)
    
    for ℓ in 1:k
        # Create submatrix with column ℓ removed
        cols_to_keep = [1:(ℓ-1); (ℓ+1):k]
        
        if length(cols_to_keep) > 0
            submatrix = B_k_x[:, cols_to_keep]
            permanents[ℓ] = BosonSampling.permanent(submatrix)
        else
            # If no columns left, permanent is 1
            permanents[ℓ] = one(ComplexF64)
        end
    end
    
    return permanents
end

#%% Cell 2: Main Algorithm Implementation

"""
    clifford_algorithm_paper(A::Matrix, n::Int)

Direct implementation of Algorithm B from the Clifford & Clifford paper.

**Input:**
- A: m×n matrix (first n columns of m×m unitary)
- n: number of photons

**Output:**
- z: sorted array of output mode indices
"""
function clifford_algorithm_paper(A::Matrix, n::Int)
    m = size(A, 1)
    
    # Line 1: r ← ∅ (EMPTY ARRAY)
    r = Int[]
    
    # Line 2: A ← PERMUTE(A) (RANDOMLY PERMUTE COLUMNS OF A)
    A = permute_columns(A, n)
    
    # Line 3: wi ← |ai,1|², i ∈ [m] (MAKE INDEXED WEIGHT ARRAY w)
    w = abs2.(A[:, 1])
    
    # Line 4: x ← SAMPLE(w) (SAMPLE INDEX x FROM w)
    x = wsample(1:m, Weights(w))
    
    # Line 5: r ← (r, x) (APPEND x TO r)
    push!(r, x)
    
    # Line 6: FOR k ← 2 TO n DO
    for k in 2:n
        # Line 7: B_k^x ← A^[k] (Submatrix of A)
        # B_k^x has rows from r and first k columns
        B_k_x = A[r, 1:k]
        
        # Line 8: COMPUTE {Per B_{k,ℓ}^x, ℓ ∈ [k]} (AS LEMMA 2)
        permanents = compute_permanents_lemma2(B_k_x, k)
        
        # Line 9: wi ← |∑_{ℓ=1}^k ai,ℓ Per B_{k,ℓ}^x|², i ∈ [m] (USING LAPLACE EXPANSION)
        w = zeros(Float64, m)
        for i in 1:m
            amplitude = zero(ComplexF64)
            for ℓ in 1:k
                amplitude += A[i, ℓ] * permanents[ℓ]
            end
            w[i] = abs2(amplitude)
        end
        
        # Line 10: x ← SAMPLE(w)
        x = wsample(1:m, Weights(w))
        
        # Line 11: r ← (r, x) (APPEND x TO r)
        push!(r, x)
    end
    # Line 12: END FOR
    
    # Line 13: z ← INCSORT(r) (SORT r IN NON-DECREASING ORDER)
    z = sort(r)
    
    # Line 14: RETURN z
    return z
end

#%% Cell 3: Interface Functions

"""
    clifford_sampler_paper(input::Input, interf::Interferometer)

Interface function for BosonSampling.jl types.
Returns occupancy vector instead of sorted mode list.
"""
function clifford_sampler_paper(input::Input, interf::Interferometer)
    # Get the unitary matrix and input arrangement
    U = interf.U
    input_arrangement = BosonSampling.fill_arrangement(input)
    
    # Extract the relevant submatrix A (first n columns)
    A = U[:, input_arrangement]
    
    # Run the algorithm
    z = clifford_algorithm_paper(A, input.n)
    
    # Convert sorted mode list to occupancy vector
    occupancy = zeros(Int, input.m)
    for mode in z
        occupancy[mode] += 1
    end
    
    return occupancy
end

#%% Cell 4: Testing and Validation

"""
    test_paper_implementation()

Test the paper implementation against our validation framework.
"""
function test_paper_implementation()
    println("="^60)
    println("TESTING CLIFFORD PAPER IMPLEMENTATION")
    println("="^60)
    
    Random.seed!(42)
    
    # Test parameters
    n, m = 2, 3
    n_samples = 5000
    
    println("System: $n photons, $m modes")
    println("Samples: $n_samples")
    
    # Create input and interferometer
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    # Generate samples using paper implementation
    println("\nGenerating samples using paper algorithm...")
    samples = []
    for i in 1:n_samples
        sample = clifford_sampler_paper(input_state, interf)
        push!(samples, sample)
        
        if i % (n_samples ÷ 10) == 0
            print(".")
        end
    end
    println(" done!")
    
    # Convert to frequency distribution
    sample_counts = Dict{Vector{Int}, Int}()
    for sample in samples
        sample_counts[sample] = get(sample_counts, sample, 0) + 1
    end
    
    # Show results
    println("\nSample distribution:")
    sorted_samples = sort(collect(sample_counts), by=x->x[2], rev=true)
    
    for (config, count) in sorted_samples
        freq = count / n_samples
        println("  $config: $(round(freq, digits=4)) ($(count) times)")
    end
    
    # Basic validation checks
    println("\nValidation checks:")
    
    # Check photon conservation
    photon_counts = [sum(sample) for sample in samples]
    conservation_ok = all(count == n for count in photon_counts)
    println("✓ Photon conservation: $(conservation_ok ? "PASS" : "FAIL")")
    
    # Check all samples are valid
    valid_samples = all(all(x >= 0 for x in sample) && length(sample) == m for sample in samples)
    println("✓ Valid samples: $(valid_samples ? "PASS" : "FAIL")")
    
    # Check diversity
    unique_count = length(unique(samples))
    println("✓ Unique configurations: $unique_count / $n_samples")
    
    println("="^60)
    
    return samples, sample_counts
end

#%% Cell 5: Run Tests

println("Algorithm B from Clifford & Clifford paper implemented!")
println("Run test_paper_implementation() to validate")

# Uncomment to run test:
# test_paper_implementation()