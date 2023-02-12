# there seems to have been strange bugs with the Bayesian validation, working sometimes only for highly symmetric interferometers but not for RandHaar
# here we try to check the validity of the samplers in both cases

n_events = 1000
n = 2
m = 2
interf = Fourier(m)
TIn = Bosonic
input_state = Input{TIn}(first_modes(n,m))

events = generate_experimental_data(n_events = n_events, n = n, m = m, interf = interf, TIn = TIn)

"""

    tvd_sampled_versus_exact_distribution(events)

From a list of events, compute the total variation distance between the sampled distribution and the exact distribution.

Returns

- the total variation distance
- the sampled distribution (normalized)
- the exact distribution
"""
function tvd_sampled_versus_exact_distribution(events::Vector{Event})

    # need to write a function that takes events and converts them into an approximate full distribution

    @argcheck length(events) > 0 "events must be non-empty"

    @argcheck assert_compatible_events(events)

    mc = MultipleCounts()
    initialise_to_empty_vectors!(mc, Int, typeof(events[1].output_measurement.s))

    for ev in events

        push!(mc.proba, 1)
        push!(mc.counts, ev.output_measurement.s)

    end

    sum_duplicates!(mc)

    to_proba!(mc)

    mc_sampled = mc

    # generating the full distribution

    ev = events[1]
    i = ev.input_state
    interf = ev.interferometer

    if typeof(ev.output_measurement.s) == ModeOccupation

        o = BosonSamplingDistribution()

    elseif typeof(ev.output_measurement.s) == ThresholdModeOccupation

        o = BosonSamplingThresholdDistribution()

    else

        error("$(typeof(ev.output_measurement.s)) not implemented")

    end

    ev_full_distribution_threshold = Event(i, o, interf)
    compute_probability!(ev_full_distribution_threshold)

    mc_full_proba = ev_full_distribution_threshold.proba_params.probability

    # convert the sampled case with the events that were not observed

    mc_sampled_all_events = deepcopy(mc_full_proba)

    for i in 1:length(mc_sampled_all_events.proba) # empty probas
        mc_sampled_all_events.proba[i] = 0
    end

    for (i, count) in enumerate(mc_sampled_all_events.counts)
        for (proba_sampled, count_sampled) in zip(mc_sampled.proba, mc_sampled.counts)
            if count == count_sampled
                mc_sampled_all_events.proba[i] = proba_sampled
            end
        end
    end

    mc_sampled_all_events

    sort!(mc_full_proba)
    sort!(mc_sampled_all_events)

    tvd(mc_full_proba, mc_sampled_all_events), mc_sampled_all_events, mc_full_proba

end

tvd_sampled_versus_exact_distribution(events)

