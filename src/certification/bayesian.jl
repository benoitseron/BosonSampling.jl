### first, a simple bayesian estimator ###

# for the theory, see 1904.12318 page 3

confidence(χ) = χ/(1+χ)

function update_confidence(event, p_q, p_a, χ)

    χ *= p_q(event)/p_a(event)
    χ

end


"""
    compute_confidence(events,p_q, p_a)

A bayesian confidence estimator: return the probability that the null hypothesis
Q is right compared to the alternative hypothesis A.
"""
function compute_confidence(events,p_q, p_a)

    confidence(compute_χ(events,p_q, p_a))
end

"""
    compute_confidence_array(events, p_q, p_a)

Return an array of the probabilities of H being true as we process more and
more events.
"""
function compute_confidence_array(events, p_q, p_a)

    χ_array = [1.]

    for event in events
        push!(χ_array, update_confidence(event, p_q, p_a, χ_array[end]))
    end

    confidence.(χ_array)

end

function compute_χ(events, p_q, p_a)
    χ = 1.
    for event in events
        χ = update_confidence(event, p_q, p_a, χ)
    end
    χ
end

"""
    compute_probability!(b::Bayesian)

Updates all probabilities associated with a `Bayesian` `Certifier`.
"""
function compute_probability!(b::Bayesian)

    b.probabilities = compute_confidence_array(b.events, b.null_hypothesis.f, b.alternative_hypothesis.f)
    b.confidence = b.probabilities[end]

end

"""
    p_B(event::Event{TIn, TOut}) where {TIn<:InputType, TOut <: FockDetection}

Outputs the probability that a given `FockDetection` would have if the `InputType` was `Bosonic` for this event.
"""
function p_B(event::Event{TIn, TOut}) where {TIn<:InputType, TOut <: FockDetection}

    interf = event.interferometer
    r = event.input_state.r
    input_state = Input{Bosonic}(r)
    output_state = event.output_measurement

    event_H = Event(input_state, output_state, interf)
    compute_probability!(event_H)

    event_H.proba_params.probability

end

"""
    p_D(event::Event{TIn, TOut}) where {TIn<:InputType, TOut <: FockDetection}

Outputs the probability that a given `FockDetection` would have if the `InputType` was `Distinguishable` for this event.
"""
function p_D(event::Event{TIn, TOut}) where {TIn<:InputType, TOut <: FockDetection}

    interf = event.interferometer
    r = event.input_state.r
    input_state = Input{Distinguishable}(r)
    output_state = event.output_measurement

    event_A = Event(input_state, output_state, interf)
    compute_probability!(event_A)

    event_A.proba_params.probability

end
