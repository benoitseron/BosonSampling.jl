"""

    full_distribution(i::Input, interf::Interferometer)
	full_distribution(params::SamplingParameters)

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

full_distribution(params::SamplingParameters) = full_distribution(params.i, params.interf)

function full_threshold_distribution(i::Input, interf::Interferometer)

    outputs = all_threshold_mode_occupations(i.n,i.m, only_photon_number_conserving = true)
    probas = zeros(length(outputs))

    @showprogress for (j, o) in enumerate(outputs)

        ev = Event(i, ThresholdFockDetection(o), interf)

        @show compute_probability!(ev)

        probas[j] = compute_probability!(ev)

    end

    if !is_lossy(interf)
        @argcheck sum(probas) ≈ 1 "The sum of the probabilities is not 1."
    end

    MultipleCounts(outputs, probas)

end

full_threshold_distribution(params::SamplingParameters) = full_threshold_distribution(params.i, params.interf)

"""

	compute_probability!(ev::Event{TIn,TOut}) where {TIn<:InputType, TOut<:BosonSamplingDistribution}

Computes the full Boson Sampling probability distribution for this event.
"""
function compute_probability!(ev::Event{TIn,TOut}) where {TIn<:InputType, TOut<:BosonSamplingDistribution}

	check_probability_empty(ev)

	ev.proba_params.precision = eps()
	ev.proba_params.failure_probability = 0

    	# ev.output_measurement.mc = full_distribution(ev.input_state, ev.interferometer)

	ev.proba_params.probability = full_distribution(ev.input_state, ev.interferometer)

end

"""

	compute_probability!(ev::Event{TIn,TOut}) where {TIn<:InputType, TOut<:BosonSamplingThresholdDistribution}

Computes the full Boson Sampling probability distribution for this event with threshold detection.
"""
function compute_probability!(ev::Event{TIn,TOut}) where {TIn<:InputType, TOut<:BosonSamplingThresholdDistribution}

	check_probability_empty(ev)

	ev.proba_params.precision = eps()
	ev.proba_params.failure_probability = 0

	# ev.output_measurement.mc = full_threshold_distribution(ev.input_state, ev.interferometer)

    ev.proba_params.probability = full_threshold_distribution(ev.input_state, ev.interferometer)

end
