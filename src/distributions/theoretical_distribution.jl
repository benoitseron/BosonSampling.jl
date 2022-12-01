"""
    theoretical_distribution(;input::Input, interf::Interferometer, i=nothing)

WARNING: use full_distribution() instead, depreciated

Compute the probability distribution of all possible output configurations of
fully/partially-indistinguishable photons through a lossless interferometer.

!!! note
    - The probabilities within the distribution are indexed following the same order as `output_mode_occupation(n,m)`
    - If `i` (with default value `nothing`) is set to an integer
      `theoretical_distribution` returns the probability to find the photons in
      the i'th configuration given by `output_mode_occupation`
"""
function theoretical_distribution(;input::Input, interf::Interferometer, i=nothing)

    input_modes = input.r.state
    number_photons = input.n
    number_modes = input.m

    if get_parametric_type(input)[1] == OneParameterInterpolation
        distinguishability = input.distinguishability_param
    else
        get_parametric_type(input)[1] == Bosonic ? distinguishability = 1.0 : distinguishability = 0.0
    end

    U = interf.U
    S = input.G.S

    input_event = fill_arrangement(input_modes)
    output_events = output_mode_occupation(number_photons, number_modes)

    function compute_pi(event)

        M = U[input_event, event]
        output_modes = fill_arrangement(event)

        if distinguishability == 1
            return abs(ryser(M)).^2
        else
            W = Array{eltype(S)}(undef, (number_photons, number_photons, number_photons))
            for rr in 1:number_photons
                for ss in 1:number_photons
                    for j in 1:number_photons
                        W[rr,ss,j] = real(M[ss,j] * conj(M[rr,j]) * S[rr,ss])
                    end
                end
            end
            return ryser_tensor(W)
        end

    end

    complete_distribution() = map(e->compute_pi(e)/factorial(number_photons), output_events)

    if i == nothing
        return complete_distribution()
    else
        return compute_pi(output_events[i])
    end

end

"""

    full_distribution(i::Input, interf::Interferometer)

Generates the complete BosonSampling distribution for the `Input` `i` and the given `Interferometer`.

"""
function full_distribution(i::Input, interf::Interferometer)

    outputs = ModeOccupation.(all_mode_configurations(i.n,i.m; only_photon_number_conserving = true))
    probas = zeros(length(outputs))

    @showprogress for (j, output) in enumerate(outputs)

        o = FockDetection(output)
        ev = Event(i, o, interf)

        probas[j] = compute_probability!(ev)

    end

    MultipleCounts(outputs, probas)

end
