function theoretical_distribution(;input::Input, distinguishability::Real, interf::Interferometer, gram_matrix::GramMatrix, i=nothing)

    input_modes = input.r.state
    number_photons = input.r.n
    number_modes = input.r.m
    U = interf.U
    S = gram_matrix.S

    input_event = fill_arrangement(input_modes)
    output_events = generate_events(number_modes, number_photons)

    function compute_pi(event)

        M = U[input_event, event]
        output_modes = fill_arrangement(event)

        if distinguishability == 1 ||distinguishability == 0
            return abs(permanent_ryser(M)).^2
        else
            return multi_dim_ryser(M, S)
        end

    end

    complete_distribution() = [compute_pi(e)/factorial(number_photons) for e in output_events]

    if i == nothing
        return complete_distribution()
    else
        return compute_pi(output_events[i])

    end

end
