"""

	compute_probability!(ev::Event{TIn,TOut}) where {TIn<:InputType, TOut<:BosonSamplingDistribution}

Computes the full Boson Sampling probability distribution for this event.
"""
function compute_probability!(ev::Event{TIn,TOut}) where {TIn<:InputType, TOut<:BosonSamplingDistribution}

	check_probability_empty(ev)

	ev.proba_params.precision = eps()
	ev.proba_params.failure_probability = 0

	ev.output_measurement.mc = full_distribution(ev.input_state, ev.interferometer)

end
