"""
Final correct implementation of Clifford algorithm

Based exactly on the unoptimized version in the main codebase,
but with proper function names and clean implementation.
"""

using LinearAlgebra
using StatsBase
using Random

"""
    clifford_sampler_final(A::Matrix, n::Int)

Final correct implementation following the exact algorithm from the codebase.
"""
function clifford_sampler_final(A::Matrix, n::Int)
    m = size(A, 1)
    
    # Randomly permute the first n columns of A
    perm = randperm(n)
    A_perm = zeros(eltype(A), m, n)
    for i in 1:n
        A_perm[:, i] = A[:, perm[i]]
    end
    A = A_perm
    
    # List of sampled output modes
    r = Vector{Int}()
    
    # Sample first photon: probability ∝ |A[i,1]|²
    weights = [abs(A[i, 1])^2 for i in 1:m]
    x = wsample(1:m, weights)
    push!(r, x)
    
    # Sample remaining photons
    for k in 2:n
        # Create submatrix B_k: rows are sampled modes r, columns are 1:k
        B_k = A[r, 1:k]
        
        # For each column l in 1:k, compute permanent with column l removed
        removed_index(l, k) = [i for i in 1:k if i != l]
        perm_array = [BosonSampling.permanent(B_k[:, removed_index(l, k)]) for l in 1:k]
        
        # For each output mode i, compute probability
        weights = [abs(sum([A[i, l] * perm_array[l] for l in 1:k]))^2 for i in 1:m]
        
        # Sample next photon
        x = wsample(1:m, weights)
        push!(r, x)
    end
    
    # Convert mode list to occupancy vector
    occupancy = zeros(Int, m)
    for mode in r
        occupancy[mode] += 1
    end
    
    return occupancy
end

"""
    clifford_sampler_input_final(input::Input, interf::Interferometer)

Final Clifford sampler for BosonSampling.jl Input and Interferometer.
"""
function clifford_sampler_input_final(input::Input, interf::Interferometer)
    A = interf.U[:, BosonSampling.fill_arrangement(input)]
    return clifford_sampler_final(A, input.n)
end