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
    cliffords_sampler(;input::Input, interf::Interferometer, nthreads=1)

Sample photons according to the [`Bosonic`](@ref) case following
[Clifford & Clifford](https://arxiv.org/pdf/2005.04214.pdf) algorithm performed
(at most) in ``O(n2^m + Poly(n,m))`` time and ``O(m)`` space.
The optional parameter `nthreads` is used to deploy the algorithm on several threads.
"""
function cliffords_sampler(;input::Input, interf::Interferometer)

    error("disabled for verification purposes")

    m = input.m
    n = input.n

    A = interf.U[:,1:n]
    A = A[:, shuffle(1:end)]

    weight_array = map(a -> abs(a).^2, A[:,1])
    sample_array = [wsample(1:m, Weights(weight_array))]

    for k in 2:n
        subA = A[:,1:k]
        @inbounds permanent_matrix = reshape(subA[sample_array, :], k-1, k)
        unormalized_pmf = LaplaceExpansion(permanent_matrix, subA)
        weight_array = unormalized_pmf/sum(unormalized_pmf)
        push!(sample_array, wsample(1:m, Weights(weight_array)))
    end

    return sort(sample_array)

end


"""
    cliffords_sampler(ev::Event{TIn, TOut}; occupancy_vector = true) where {TIn<:InputType, TOut <: FockSample}

Sampler for an [`Event`](@ref). Note the difference of behaviour if `occupancy_vector = true`.
"""
function cliffords_sampler(ev::Event{TIn, TOut}; occupancy_vector = true) where {TIn<:InputType, TOut <: FockSample}

    error("disabled for verification purposes")

    s = cliffords_sampler(input = ev.input_state, interf = ev.interferometer)

    occupancy_vector ? mode_occupancy_to_occupancy_vector(s, ev.input_state.m) : s
end
