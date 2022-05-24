"""
    cliffords_sampler(;input::Input, interf::Interferometer)

Sample photons according to the [`Bosonic`](@ref) case following
[Clifford & Clifford](https://arxiv.org/pdf/2005.04214.pdf) algorithm performed
(at most) in ``O(n2^m + Poly(n,m))`` time and ``O(m)`` space.
"""
function cliffords_sampler(;input::Input, interf::Interferometer)

    m = input.m
    n = input.n
    perm_n_ = permutations(1:n)

    U = interf.U
    A = U[:, 1:n]
    A = A[:, nthperm(1:n, rand(1:factorial(big(n))))]

    weight_array = [abs(A[i,1])^2 for i in 1:m]
    sample_array = [wsample(1:m, Weights(weight_array))]

    for k in 2:n
        k_ = collect(1:k)
        permanent_matrix = reshape(A[sample_array, k_], length(sample_array), k)
        unormalized_pmf = LaplaceExpansion(permanent_matrix, A[:,k_])
        weight_array = unormalized_pmf/sum(unormalized_pmf)
        push!(sample_array, wsample(1:m, Weights(weight_array)))
    end

    return sort(sample_array)

end
