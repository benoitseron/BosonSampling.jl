"""
    gaussian_sampler(;input::GaussianInput{T}, nsamples=Int(1e3), burn_in=200, thinning_rate=100) where {T<:Gaussian}

Simulate noiseless Gaussian Boson Sampling with photon number resolving detectors via a MIS sampler.
!!! note Reference
    [The Boundary for Quantum Advantage in Gaussian Boson Sampling](https://arxiv.org/pdf/2108.01622.pdf)
"""
function gaussian_sampler(;input::GaussianInput{T}, nsamples=Int(1e3), burn_in=200, thinning_rate=100) where {T<:Gaussian}

    Ri = input.displacement
    V = input.covariance_matrix
    n = div(LinearAlgebra.checksquare(V), 2)
    D, S = williamson_decomp(V)
    id = Matrix{eltype(V)}(I, size(V))
    Did = D - id
    for i in 1:size(Did)[1]
        for j in 1:size(Did)[2]
            Did[i,j] < 0 ? Did[i,j] = 0 : nothing
        end
    end

    Tv = S * transpose(S)
    W = S * Did * transpose(S)
    sqrtW = real.(S * sqrt.(Did))

    A = A_mat(Tv)
    Q = husimiQ_matrix(Tv)
    detQ = det(Q)
    B = A[1:n, 1:n]
    B2 = abs.(B).^2

    sample_R() = sqrtW * rand(MvNormal(Ri))
    compute_displacement(R::AbstractVector) = [(R[i]+im*R[i+1])/sqrt(2) for i in 1:2:length(R)-1]
    sample_displacement() = compute_displacement(sample_R())

    function sample_pattern()

        displacement = sample_displacement()
        G = abs.(displacement).^2
        pattern = [rand(Poisson(g)) for g in G]

        for j in 1:length(displacement)
            pattern[j] += 2 * rand(Poisson(0.5 * B2[j,j]))
            for k in j+1:length(displacement)
                m = rand(Poisson(B2[j,k]))
                pattern[j] += m
                pattern[k] += m
            end
        end

        return pattern, displacement

    end

    function valid_sampling(N=nothing)

        if N == nothing
            return sample_pattern()
        else
            pattern = []
            displacement = []
            count = nothing

            while count != N
                pattern, displacement = sample_pattern()
                count = sum(pattern)
            end

            return pattern, displacement
        end

    end

    function sampled_probability(pattern, displacement)

        G = abs.(displacement).^2
        Gn = duplicate_row_col(G, pattern)

        Bn = duplicate_row_col(B2, pattern)
        for i in 1:length(G)
            Bn[i,i] = Gn[i]
        end

        lhaf = hafnian(Bn, loop=true)
        factor = exp(-sum(G))
        rescale = prod(factorial.(pattern))

        return real((factor*lhaf)/rescale)

    end

    function target_probability(pattern, displacement)

        gamma = displacement .- B*conj(displacement)
        gamma_n = duplicate_row_col(gamma, pattern)

        Bn = duplicate_row_col(B, pattern)
        for i in 1:length(gamma_n)
            Bn[i,i] = gamma_n[i]
        end

        lhaf = hafnian(Bn, loop=true)
        arg_exp = -1/2 * (norm(displacement)^2 - dot(conj(displacement), B*conj(displacement)))
        factor = abs(exp(arg_exp))^2
        rescale = prod(factorial.(pattern)) * sqrt(detQ)

        return real(factor * lhaf / rescale)

    end

    pattern, displacement = valid_sampling(4)
    patterns = [pattern]
    pattern_chain = pattern
    p_proposal_chain = sampled_probability(pattern, displacement)
    p_target_chain = target_probability(pattern, displacement)

    pattern_chain_ = []
    outcomes = []
    p_proposal_ = []
    p_target_ = []

    for i in 1:nsamples

        pattern, displacement = valid_sampling(4)
        p_proposal = sampled_probability(pattern, displacement)
        p_target = target_probability(pattern, displacement)
        p_accept = min(1, (p_proposal_chain*p_target)/(p_target_chain*p_proposal))

        if rand() <= p_accept
            push!(pattern_chain_, pattern)
            p_proposal_chain = p_proposal
            p_target_chain = p_target
            pattern_chain = pattern
            push!(outcomes, 1)
        else
            push!(pattern_chain_, pattern_chain)
            push!(outcomes, 0)
        end

        push!(p_proposal_, p_proposal)
        push!(p_target_, p_target)

    end

    res = []
    for i in thinning_rate:nsamples
        if outcomes[i] == 1
            push!(res, pattern_chain_[i])
        end
    end

    return res

end
