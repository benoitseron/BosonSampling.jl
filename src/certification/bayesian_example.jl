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

n_events = 500
n = 4
m = 128
interf = Hadamard(m)
TIn = Bosonic
input_state = Input{TIn}(first_modes(n,m))

    # n = 2
    # sparsity = 3
    # m = sparsity * n

    # # x = 0.9
    # # T = OneParameterInterpolation
    # TIn = Bosonic
    # mode_occ = equilibrated_input(sparsity, m)

    # d = Uniform(0,2pi)
    # ϕ = nothing # rand(d,m)
    # η_loss_lines =  0.86 * ones(m)
    # η_loss_bs = 0.93 * ones(m-1)
    # η_loss_source = get_η_loss_source(m, QuantumDot(13.5/80))

    # η = rand(m-1)

    # params = LoopSamplingParameters(n=n, m=m, η = η, η_loss_bs = η_loss_bs, η_loss_lines = η_loss_lines, η_loss_source = η_loss_source, ϕ = ϕ,T=T, mode_occ = mode_occ)

    # interf = build_loop!(params)
    # input_state = Input{TIn}(first_modes(n,m))

    


interf

events = []

for i in 1:n_events

    # note that we don't compute the event probability
    # as we would just have experimental observations
    # of counts

    ev = Event(input_state, FockSample(), interf)
    BosonSampling.sample!(ev)

    ev = convert(Event{TIn, FockDetection}, ev)

    push!(events, ev)

end


# now we have the vector of observed events with probabilities

events

events[100].interferometer

# next, from events, recover the probabilities under both
# hypothesis for instance

p_B(events[1])
p_D(events[1])

# hypothesis : the events were from a bosonic distribution
# for that we use the Bayesian type

p_q = HypothesisFunction(p_B)
p_a = HypothesisFunction(p_D)

certif = Bayesian(events, p_q, p_a)
BosonSampling.certify!(certif, max_χ = 10000)
certif.confidence

scatter(certif.probabilities)

######## it looks like everything goes bad as soon as using RandHaar while it works for Fourier, Hadamard... it looks like the sampler is bad... maybe not even sampling but having random_mode_occupation_collisionless or something alike
####### this doesn't seem to impact the partition validator nor the full bunching



###### Bayesian tests for partitions ######

# need to provide the interferometer as well as some data to compute the probabilities
# this can be extracted from experimental if given it but here we have to generate it

# generate what would be experimental data
m = 14
n = 2

TIn = Bosonic
mode_occ = equilibrated_input(sparsity, m)

d = Uniform(0,2pi)
ϕ = nothing # rand(d,m)
η_loss_lines =  nothing#0.86 * ones(m)
η_loss_bs = nothing#0.93 * ones(m-1)
η_loss_source = nothing#get_η_loss_source(m, QuantumDot(13.5/80))

η = rand(m-1)

params = LoopSamplingParameters(n=n, m=m, η = η, η_loss_bs = η_loss_bs, η_loss_lines = η_loss_lines, η_loss_source = η_loss_source, ϕ = ϕ,T=T, mode_occ = mode_occ)

interf = build_loop!(params)



events = generate_experimental_data(n_events = 1000, n = n,m = m, interf = build_loop!(params), TIn = Bosonic)


n_subsets = 4
part = equilibrated_partition(m,n_subsets)
certif = BayesianPartition(events, Bosonic(), Distinguishable(), is_lossy(interf) ? to_lossy(part) : part)

certify!(certif, max_χ = Inf, min_χ  = 0.0001)
plot(certif.probabilities)


# full bunching

subset_size = m-n
subset = Subset(first_modes(subset_size, m))

fb = FullBunching(events, Bosonic(), Distinguishable(), subset_size)

certify!(fb)
