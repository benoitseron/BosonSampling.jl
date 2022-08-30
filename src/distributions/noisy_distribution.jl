"""
    noisy_distribution(;input::Input, loss::Real, interf::Interferometer, exact=true, approx=true, samp=true; error=1e-4, failure_probability=1e-4)

Compute the exact and/or approximated and/or sampled probability distribution of
all possible output configurations of partially-distinguishable photons through a
lossy interferometer.
By default, `exact`, `approx` and `samp` are set to `true` meaning that `noisy_distribution`
returns an array containing the three distributions.

!!! note
    - The probabilities within a distribution are indexed following the same order as [`output_mode_occupation(n,m)`](@ref)
    - The approximated distribution has error and failure probability of ``1e^{-4}``.
!!! note "Reference"
    [https://arxiv.org/pdf/1809.01953.pdf](https://arxiv.org/pdf/1809.01953.pdf)
"""
function noisy_distribution(;input::Input, loss::Real, interf::Interferometer, exact=true, approx=true, samp=true; error=1e-4, failure_probability=1e-4)

    output = []
    ϵ = error
    δ = failure_probability

    input_modes = input.r.state
    number_photons = input.n
    number_modes = input.m
    U = interf.U

    if get_parametric_type(input)[1] == OneParameterInterpolation
        distinguishability = input.distinguishability_param
    else
        get_parametric_type(input)[1] == Bosonic ? distinguishability = 1.0 : distinguishability = 0.0
    end

    number_output_photons = trunc(Int, number_photons*loss)
    input_occupancy_modes = fill_arrangement(input_modes)

    if loss == 1 ||distinguishability == 0
        throw(ArgumentError("invalid input parameters"))
    end

    ki = (log(ϵ*δ*(1-loss*distinguishability^2)/2))/log(loss*distinguishability^2)
    ki = 10^(trunc(Int, log(10, ki)))
    ki < 10 ? k = number_photons : k = min(ki, number_photons)
    kmax = k-1

    combs = collect(combinations(input_occupancy_modes, k))
    nlist = output_mode_occupation(k, number_modes)

    function compute_pi(i::Integer, ans)

        res = 0
        for j = 1:length(combs)
            perm = collect(permutations(combs[j]))
            for l = 1:length(perm)

                count = sum(collect(Int(perm[l][ll] == combs[j][ll]) for ll = 1:length(perm[l])))
                korder = k - count

                if ans == true
                    pterm = distinguishability^korder * ryser(U[combs[j], nlist[i]] .* conj(U[perm[l], nlist[i]]))
                    pterm = real(pterm / (binomial(number_photons, k) * factorial(k)))
                    res += pterm
                elseif ans == false && korder <= kmax
                    pterm = distinguishability^korder * ryser(U[combs[j], nlist[i]] .* conj(U[perm[l], nlist[i]]))
                    pterm = real(pterm / (binomial(number_photons, k) * factorial(k)))
                    res += pterm
                end

            end
        end

        return res

    end

    compute_full_distribution(ans) = map(i->compute_pi(i,ans), 1:length(nlist))

    function sampling(n_samples = 1e5)

        b = nlist[1]
        comb = input_occupancy_modes[1:k]
        proba = compute_pi(1, true)
        current_pos = 1
        lsample = []
        bprob = zeros(length(nlist))

        for i = 1:n_samples

            pos_test = rand(1:length(nlist))
            b_test = nlist[pos_test]
            comb_test = combs[rand(1:length(combs))]
            perm = collect(permutations(comb_test))
            prob_test = 0

            for j = 1:length(perm)
                count = sum(collect(Int(k - sum(perm[j]) == comb_test)))
                if count <= kmax
                    prob_test += ryser(U[comb_test, b_test] .* conj(U[perm[j], b_test]))
                end
            end

            prob_test = real(prob_test)

            if prob_test > proba
                b = b_test
                proba = prob_test
                comb = comb_test
                current_pos = pos_test
            elseif prob_test / proba > rand()
                b = b_test
                proba = prob_test
                comb = comb_test
                current_pos = pos_test
            end

            bprob[current_pos] += 1 / n_samples
            push!(lsample) = current_pos

        end

        return bprob

    end

    exact ? push!(output, compute_full_distribution(true)) : nothing
    approx ? push!(output, compute_full_distribution(false)) : nothing
    samp ? push!(output, sampling()) : nothing

    return output

end
