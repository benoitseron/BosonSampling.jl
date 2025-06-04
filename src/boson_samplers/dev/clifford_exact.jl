"""
Exact implementation of the Clifford algorithm following the codebase precisely
"""

using LinearAlgebra
using StatsBase
using Random

"""
    clifford_exact(A, n; occupancy_vector = true)

Exact copy of the unoptimized Clifford algorithm from the main codebase.
"""
function clifford_exact(A, n; occupancy_vector = true)
    # most basic sampler
    m = size(A,1)
    
    r = Vector{Int}()
    
    # randomly permute the first n columns of A
    function permute_columns(A,n)
        perm = randperm(n)
        B = zeros(eltype(A),m,n)
        for i in 1:n
            B[:,i] = A[:,perm[i]]
        end
        B
    end
    
    A = permute_columns(A,n)
    
    weights = [abs(A[i,1])^2 for i in 1:m]
    
    x = wsample(1:m, weights)
    
    push!(r,x)
    
    indexes_remove(r,k) = [i for i in 1:k if i ∉ r]
    
    for k in 2:n
        B_k = A[r, 1:k]
        
        removed_index(l,k) = [i for i in 1:k if i != l]
        perm_array = [BosonSampling.permanent(B_k[:,removed_index(l,k)]) for l in 1:k]
        
        weights = [abs(sum([A[i,l] * perm_array[l] for l in 1:k]))^2 for i in 1:m]
        
        x = wsample(1:m, weights)
        
        push!(r,x)
    end
    
    occupancy_vector ? BosonSampling.mode_occupancy_to_occupancy_vector(r, m) : r
end

"""
    clifford_sampler_exact(input::Input, interf::Interferometer)

Exact wrapper using the same interface as the original.
"""
function clifford_sampler_exact(input::Input, interf::Interferometer)
    A = interf.U[:, BosonSampling.fill_arrangement(input)]
    return clifford_exact(A, input.n, occupancy_vector = true)
end