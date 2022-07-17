### tests with standard boson sampling ###

# we will test that we are indeed in the bosonic case
# compared to the distinguishable one

# first we generate a series of bosonic events

n_events = 20
n = 3
m = 8
interf = RandHaar(m)
input_state = Input{Bosonic}(first_modes(n,m))

events = []

for i in 1:n_events

    # generate a random output pattern
    output_state = FockDetection(random_mode_occupation(n,m))

    # note that we don't compute the event probability
    # as we would just have experimental observations
    # of counts

    this_event = Event(input_state, output_state, interf)
    push!(events, this_event)

end

# now we have the vector of observed events with probabilities

events

# next, from events, recover the probabilities under both
# hypothesis for instance

p_B(events[1])
p_D(events[1])

# hypothesis : the events were from a bosonic distribution
# for that we use the Bayesian type

p_q = HypothesisFunction(p_B)
p_a = HypothesisFunction(p_D)

certif = Bayesian(events, p_q, p_a)
compute_probability!(certif)
certif.confidence

scatter(certif.probabilities)
