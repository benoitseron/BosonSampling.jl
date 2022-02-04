function cliffords_sampler(;input::Input, interf::Interferometer)

    m = input.r.m
    n = input.r.n
    n_ = collect(1:n)
    perm_n = collect(permutations(n_))

    A = interf.U[:,n_]
    A = A[:,rand(perm_n)]

    weight_array = [abs(A[i,1])^2 for i in 1:m]
    sample_array = [wsample(1:m, Weights(weight_array))]

    for k in 2:n
        k_ = collect(1:k)
        permanent_matrix = A[sample_array,k_]
        unormalized_weight_array = LaplaceExpansion(permanent_matrix, A[:,k_])
        weight_array = unormalized_weight_array/sum(unormalized_weight_array)
        push!(sample_array, wsample(1:m, Weights(weight_array)))
    end

    return sort(sample_array)

end
