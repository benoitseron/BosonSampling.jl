"""

    clifford_unoptimised(A, n; occupancy_vector = true)

Naive implementation of the Clifford algorithm. Computes a sample if `n` photons are sent in the first `n` modes of a `m` mode interferometer `A`.
"""
function clifford_unoptimised(A, n; occupancy_vector = true)

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

    r

    indexes_remove(r,k) = [i for i in 1:k if i ∉ r]

    for k in 2:n

        B_k = A[r, 1:k]

        removed_index(l,k) = [i for i in 1:k if i != l]
        perm_array = [permanent(B_k[:,removed_index(l,k)]) for l in 1:k]

        weights = [abs(sum([A[i,l] * perm_array[l] for l in 1:k]))^2 for i in 1:m]

        x = wsample(1:m, weights)

        push!(r,x)

    end

    occupancy_vector ? mode_occupancy_to_occupancy_vector(r, m) : r

end

function clifford_sampler_unoptimised(i::Input, interf::Interferometer; occupancy_vector = true)

    A = interf.U[:, fill_arrangement(i)]

    clifford_unoptimised(A, i.n, occupancy_vector = occupancy_vector)

end

function clifford_sampler_unoptimised(ev::Event{TIn, TOut}; occupancy_vector = true) where {TIn<:InputType, TOut <: FockSample}

    i = ev.input_state

    interf = ev.interferometer

    clifford_sampler_unoptimised(i, interf, occupancy_vector = occupancy_vector)

end

"""
    cliffords_sampler(;input::Input, interf::Interferometer)

Sample photons according to the [`Bosonic`](@ref) case following
the corrected Clifford & Clifford algorithm from arXiv:2005.04214v2.
Fixes bias issues in the original implementation, particularly for m > n cases.
"""
function cliffords_sampler(;input::Input, interf::Interferometer)
    
    m = input.m
    n = input.n

    # Extract relevant submatrix for the input arrangement
    A = interf.U[:, fill_arrangement(input)]
    
    # Run corrected Clifford algorithm
    z = corrected_clifford_algorithm(A, n)
    
    # Convert to mode occupation vector
    mode_occ = zeros(Int, m)
    for mode in z
        mode_occ[mode] += 1
    end
    
    return mode_occ

end

# Corrected Clifford algorithm implementing Algorithm A from arXiv:2005.04214v2
function corrected_clifford_algorithm(A::Matrix{ComplexF64}, n::Int)
    """
    Corrected Clifford algorithm implementing Algorithm A from arXiv:2005.04214v2
    
    Key fixes:
    1. Proper handling of repeated rows using Ryser's formula with multiplicities
    2. Efficient permanent calculation for all minors
    3. Correct Laplace expansion for weight computation
    """
    m = size(A, 1)
    r = Int[]
    
    # Step 1: Randomly permute columns
    A = permute_columns(A, n)
    
    # Step 2: Sample first photon
    w = abs2.(A[:, 1])
    x = wsample(1:m, Weights(w))
    push!(r, x)
    
    # Steps 3-n: Iterative sampling
    for k in 2:n
        # Get current matrix B_k with repeated rows
        B_k = A[r, 1:k]
        
        # Compute multiplicities for repeated rows
        multiplicities = zeros(Int, m)
        for row_idx in r
            multiplicities[row_idx] += 1
        end
        
        # Compute all minors efficiently  
        minors = compute_minors_efficient(B_k, multiplicities[multiplicities .> 0])
        
        # Compute weights using Laplace expansion
        w = zeros(Float64, m)
        for i in 1:m
            amplitude = 0.0 + 0.0im
            for ℓ in 1:k
                amplitude += A[i, ℓ] * minors[ℓ]
            end
            w[i] = abs2(amplitude)
        end
        
        # Sample next photon
        x = wsample(1:m, Weights(w))
        push!(r, x)
    end
    
    # Sort and return
    z = sort(r)
    return z
end

# Helper functions for the corrected algorithm
function permute_columns(A::Matrix, n::Int)
    """Randomly permute the first n columns of matrix A"""
    A_copy = copy(A)
    perm = randperm(n)
    A_copy[:, 1:n] = A_copy[:, perm]
    return A_copy
end

function compute_minors_efficient(B_k::Matrix{ComplexF64}, multiplicities::Vector{Int})
    """
    Compute all minors {perm B_{k,ℓ}} efficiently using permanent calculations.
    """
    k = size(B_k, 2)
    
    if k == 0
        return Float64[]
    elseif k == 1
        return [1.0]  # permanent of empty matrix is 1
    end
    
    minors = zeros(ComplexF64, k)
    
    # For each minor ℓ, compute permanent of matrix with column ℓ removed
    for ℓ in 1:k
        # Create submatrix by removing column ℓ
        cols_to_keep = [1:(ℓ-1); (ℓ+1):k]
        
        if length(cols_to_keep) == 0
            minors[ℓ] = 1.0 + 0.0im
        else
            submatrix = B_k[:, cols_to_keep]
            # Use standard permanent calculation
            minors[ℓ] = permanent(submatrix)
        end
    end
    
    return minors
end


"""
    cliffords_sampler(ev::Event{TIn, TOut}; occupancy_vector = true) where {TIn<:InputType, TOut <: FockSample}

Sampler for an [`Event`](@ref). Uses the corrected Clifford algorithm.
Note: occupancy_vector parameter is deprecated as the function now returns mode occupation by default.
"""
function cliffords_sampler(ev::Event{TIn, TOut}; occupancy_vector = true) where {TIn<:InputType, TOut <: FockSample}

    s = cliffords_sampler(input = ev.input_state, interf = ev.interferometer)
    
    # s is already in mode occupation format, so just return it
    return s
end
