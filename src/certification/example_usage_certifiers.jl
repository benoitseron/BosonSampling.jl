using Revise

using BosonSampling:sample!
using BosonSampling
using Plots
using ProgressMeter
using Distributions
using Random
using Test
using ArgCheck
using StatsBase
using ColorSchemes
using Interpolations
using Dierckx
using LinearAlgebra
using PrettyTables
using LaTeXStrings
using JLD
using AutoHashEquals
using HypothesisTests


### tests with standard boson sampling ###

# we will test that we are indeed in the bosonic case
# compared to the distinguishable one

# first we generate a series of bosonic events

n = 8
m = 8
interf = RandHaar(m)
TIn = Bosonic

events = generate_experimental_data(n_events = 100, n = n,m = m, interf = RandHaar(m), TIn = Bosonic)


# next, from events, recover the probabilities under both
# hypothesis for instance

p_B(events[1])
p_D(events[1])

# hypothesis : the events were from a bosonic distribution
# for that we use the Bayesian type

p_q = HypothesisFunction(p_B)
p_a = HypothesisFunction(p_D)

certif = Bayesian(events, p_q, p_a)
BosonSampling.certify!(certif)
certif.confidence

scatter(certif.probabilities)

###### Bayesian tests for partitions ######

# need to provide the interferometer as well as some data to compute the probabilities
# this can be extracted from experimental if given it but here we have to generate it

# generate what would be experimental data
m = 14
n = 5
events = generate_experimental_data(n_events = 10, n = n,m = m, interf = RandHaar(m), TIn = Bosonic)


n_subsets = 3
part = equilibrated_partition(m,n_subsets)
certif = BayesianPartition(events, Bosonic(), Distinguishable(), part)

certify!(certif)
plot(certif.probabilities)

part

### full bunching###

m = 20
n = 5
events = generate_experimental_data(n_events = 1000, n = n,m = m, interf = RandHaar(m), TIn = Bosonic)

subset_size = m-n
subset = Subset(first_modes(subset_size, m))

fb = FullBunching(events, Bosonic(), Distinguishable(), subset_size)

certify!(fb)

### correlators ###

m = 20
n = 5
interf = RandHaar(m)
events = generate_experimental_data(n_events = 1000, n = n,m = m, interf = interf, TIn = Bosonic)

correlators_nm_cv_s(interf, input_state)
