function noisy_distribution(;input, distinguishability, reflectivity, interf, exact=true, approx=true, samp=true)

    # https://arxiv.org/pdf/1809.01953.pdf

    @warn "changed ryser->naive since actual ryser does not support non square matrix"

    output = []
    ϵ = 1e-4
    δ = 1e-4

    input_modes = input.r.state
    number_photons = input.r.n
    number_modes = input.r.m
    U = interf.U

    number_output_photons = trunc(Int, number_photons*reflectivity)
    input_occupancy_modes = fill_arrangement(input_modes)

    if reflectivity == 1 ||distinguishability == 0
        throw(ArgumentError("invalid input parameters"))
    end

    ki = (log(ϵ*δ*(1-reflectivity*distinguishability^2)/2))/log(reflectivity*distinguishability^2)
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

                count = sum(
                    collect(
                        Int(perm[l][ll] == combs[j][ll]) for
                        ll = 1:length(perm[l])
                    ),
                )
                korder = k - count

                if ans == true

                    pterm =
                        distinguishability^korder * ryser(
                            U[combs[j], nlist[i]] .* conj(U[perm[l], nlist[i]]),
                        )
                    pterm = real(
                        pterm / (binomial(number_photons, k) * factorial(k)),
                    )
                    res += pterm

                elseif ans == false && korder <= kmax

                    pterm =
                        distinguishability^korder * ryser(
                            U[combs[j], nlist[i]] .* conj(U[perm[l], nlist[i]]),
                        )
                    pterm = real(
                        pterm / (binomial(number_photons, k) * factorial(k)),
                    )
                    res += pterm
                end
            end
        end

        return res

    end

    compute_full_distribution(ans) = [compute_pi(i,ans) for i = 1:length(nlist)]

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
                    prob_test += ryser(
                        U[comb_test, b_test] .* conj(U[perm[j], b_test]),
                    )
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
