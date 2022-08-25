# generate what would be experimental data
events = generate_experimental_data(n_events = 400, n = 3,m = 10, interf = RandHaar(10), TIn = Bosonic)

# hypothesis : the events were from a bosonic distribution

p_q = HypothesisFunction(p_B)
p_a = HypothesisFunction(p_D)

# use Bayesian certification
certif = Bayesian(events, p_q, p_a)
BosonSampling.certify!(certif)
certif.confidence

scatter(certif.probabilities)
