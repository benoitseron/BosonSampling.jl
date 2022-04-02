### first, a simple bayesian estimator ###

confidence(χ) = χ/(1+χ)

function update_confidence(event, p_q, p_a, χ)

    χ *= p_q(event)/p_a(event)
    χ

end

function compute_χ(events, p_q, p_a)
    χ = 1.
    for event in events
        χ = update_confidence(event, p_q, p_a, χ)
    end
    χ
end

function compute_confidence(events,p_q, p_a)

    """a bayesian confidence estimator:

    returns the probability that the null hypothesis Q is right compared to
    the alternative hypothesis A

    the functions take in events
    p_q is a function that takes in an event and gives its probability in the null hypothesis
    p_a is the same for the alternative hypothesis
    χ is the ratio between null hypothesis probabilities and alternative ones

    """

    confidence(compute_χ(events,p_q, p_a))
end

function compute_confidence_array(events, p_q, p_a)

    """gives an array of the probability of H being true
    as we process more and more events"""

    χ_array = [1.]

    for event in events
        push!(χ_array, update_confidence(event, p_q, p_a, χ_array[end]))
    end

    confidence.(χ_array)

end

### tests with standard boson sampling ###

# we will test that we are indeed in the bosonic case
# compared to the distinguishable one

# first we generate a series of bosonic events

n_events = 10
n = 3
m = 4
interf = RandHaar(m)
input_state = Input{Bosonic}(first_modes(n,m))

events = []

for i in 1:n_events

    # generate a random output pattern
    output_state = OutputMeasurement{FockDetection}(random_mode_occupation_collisionless(n,m))

    # compute the event probability
    this_event = Event(input_state, output_state, interf)
    compute_probability!(this_event)
    push!(events, this_event)

end

# now we have the vector of observed events with probabilities

events

# next, from events, recover the probabilities under both
# hypothesis

function p_B(event::Event)

    interf = event.interferometer
    r = event.input_state.r
    input_state = Input{Bosonic}(r)
    output_state = event.output_measurement

    event_H = Event(input_state, output_state, interf)
    compute_probability!(event_H)

    event_H.proba_params.probability

end

function p_D(event::Event)

    interf = event.interferometer
    r = event.input_state.r
    input_state = Input{Distinguishable}(r)
    output_state = event.output_measurement

    event_A = Event(input_state, output_state, interf)
    compute_probability!(event_A)

    event_A.proba_params.probability

end

# hypothesis : the events were from a bosonic distribution

p_q = p_B
p_a = p_D

confidence(compute_χ(events, p_q, p_a))

# hypothesis : the events were from a distinguishable distribution

p_q = p_D
p_a = p_B

confidence(compute_χ(events, p_q, p_a))

@test confidence(compute_χ(events, p_q, p_a)) + confidence(compute_χ(events, p_a, p_q)) ≈ 1 atol = 1e-6

##### the only thing I see is a problem in the probabilities themselves?
