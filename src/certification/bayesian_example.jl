using Revise

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



### tests with standard boson sampling ###

# we will test that we are indeed in the bosonic case
# compared to the distinguishable one

# first we generate a series of bosonic events

n_events = 200
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

###### Bayesian tests for partitions ######

# need to provide the interferometer as well as some data to compute the probabilities
# this can be extracted from experiemental if given it but here we have to generate it

n_events = 100
n = 5
m = 14
interf = RandHaar(m)

# we generate the experimental data

input_state = Input{Distinguishable}(first_modes(n,m))

events = []

for i in 1:n_events

    # note that we don't compute the event probability
    # as we would just have experimental observations
    # of counts

    ev = Event(input_state, FockSample(), interf)
    sample!(ev)

    push!(events, ev)

end


# we now compute the partition probabilities for the hypothesis of Bosonic and Distinguishable inputs

ib = Input{Bosonic}(first_modes(n,m))
id = Input{Distinguishable}(first_modes(n,m))
n_subsets = 3

part = equilibrated_partition(m,n_subsets)
o = PartitionCountsAll(part)

evb = Event(ib,o,interf)
evd = Event(id,o,interf)

pb = compute_probability!(evb)
pd = compute_probability!(evd)

p_partition_B(ev) = p_partition(to_partition_count(ev, part), evb)
p_partition_D(ev) = p_partition(to_partition_count(ev, part), evd)

p_q = HypothesisFunction(p_partition_B)
p_a = HypothesisFunction(p_partition_D)

certif = Bayesian(events, p_q, p_a)
compute_probability!(certif)
certif.confidence

scatter(certif.probabilities)
