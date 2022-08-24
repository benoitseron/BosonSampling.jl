"""
    cliffords_sampler(;input::Input, interf::Interferometer, nthreads=1)

Sample photons according to the [`Bosonic`](@ref) case following
[Clifford & Clifford](https://arxiv.org/pdf/2005.04214.pdf) algorithm performed
(at most) in ``O(n2^m + Poly(n,m))`` time and ``O(m)`` space.
The optional parameter `nthreads` is used to deploy the algorithm on several threads.
"""
function cliffords_sampler(;input::Input, interf::Interferometer)

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
    s = cliffords_sampler(input = ev.input_state, interf = ev.interferometer)

    occupancy_vector ? mode_occupancy_to_occupancy_vector(s, ev.input_state.m) : s
end
