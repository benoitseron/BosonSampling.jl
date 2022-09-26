# generate what would be experimental data
m = 14
n = 5
events = generate_experimental_data(n_events = 1000, n = n,m = m, interf = RandHaar(m), TIn = Bosonic)

# define a certification protocol
# using binning into a 3-partition
# with null hypothesis a bosonic input
# and alternative a distinguishable input
n_subsets = 3
part = equilibrated_partition(m,n_subsets)
certif = BayesianPartition(events, Bosonic(), Distinguishable(), part)

# compute the level of confidence
certify!(certif)

# plot the bayesian confidence over sample number
plot(certif.probabilities)
