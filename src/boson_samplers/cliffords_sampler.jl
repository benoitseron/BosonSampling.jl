function cliffords_sampler(;input::Input, interf::Interferometer)

    m = input.m
    n = input.n
    n_ = collect(1:n)
    perm_n_ = collect(permutations(n_))

    U = interf.U
    A = U[:, n_]
    A = A[:, rand(perm_n_)]

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
