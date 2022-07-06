using Revise
using BosonSampling
using Permanents
using PrettyTables
using ArgCheck
using Test
using Permanents
using Plots
using Test
using Combinatorics
using Random
using IterTools
using Statistics
using LinearAlgebra #removed so as to be able to use generic types such as BigFloats, can be put back if needed
using PolynomialRoots
using StatsBase
#using JLD
using CSV
using DataFrames
using Tables
using Plots#; plotly() #plotly allows to make "dynamic" plots where you can point the mouse and see the values, they can also be saved like that but note that this makes them heavy (and html instead of png)
using PrettyTables #to display large matrices
using Roots
using BenchmarkTools
using Optim
using ProgressMeter
using Parameters
using ArgCheck
using Distributions
using Luxor
using UnPack
# using Distributed
# using SharedArrays

const ATOL = 1e-10
#
# using DocumenterTools

### first, a simple bayesian estimator ###

"""
    confidence(χ)

Fractional certainty associated with a bayesian ratio χ.
"""
confidence(χ) = χ/(1+χ)


"""
    update_confidence(event, p_q, p_a, χ)

Updates the bayesian ratio χ given an event, the probabilities of the two
hypothses and the previous χ.
"""
function update_confidence(event, p_q, p_a, χ)

    χ *= p_q(event)/p_a(event)
    χ

end

"""
    compute_χ(event, p_q, p_a, χ)

Computes the bayesian ratio χ given a list of events, and how to compute the probabilities of the two hypotheses for each event of the list of events.
"""
function compute_χ(events, p_q, p_a)
    χ = 1.
    for event in events
        χ = update_confidence(event, p_q, p_a, χ)
    end
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

### tests with standard boson sampling ###

# we will test that we are indeed in the bosonic case
# compared to the distinguishable one

# first we generate a series of bosonic events

n_events = 10
n = 8
m = 10
interf = RandHaar(m)
input_state_bos = Input{Bosonic}(first_modes(n,m))
input_state_dist = Input{Distinguishable}(first_modes(n,m))

input_states = [input_state_bos, input_state_dist]
######### changing to Input{Distinguishable}(first_modes(n,m)) doesn't change the output result although it very much should!

events_bos = []
events_dist = []

events_list = (events_bos, events_dist)

for i in 1:n_events

    # generate a random output pattern
    output_state = FockDetection(random_mode_occupation_collisionless(n,m))

    # compute the event probability

    for (input_state, events) in zip(input_states, events_list)
        this_event = Event(input_state, output_state, interf)
        compute_probability!(this_event)
        push!(events, this_event)
    end


end

# now we have the vector of observed events with probabilities

events_bos


# next, from events, recover the probabilities under both
# hypothesis

"""
    p_B(event::Event)

Computes the probability associated a given even configuration (input and outpute state mode occupancy with a given `Interferometer`) if the input was `Bosonic`.
"""
function p_B(event::Event)

    interf = event.interferometer
    r = event.input_state.r
    input_state = Input{Bosonic}(r)
    output_state = event.output_measurement

    event_H = Event(input_state, output_state, interf)
    compute_probability!(event_H)

    event_H.proba_params.probability

end

"""
    p_D(event::Event)

Computes the probability associated a given even configuration (input and outpute state mode occupancy with a given `Interferometer`) if the input was `Distinguishable`.
"""
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

p_q

confidence(compute_χ(events, p_q, p_a))


# hypothesis : the events were from a distinguishable distribution

p_q = p_D
p_a = p_B

p_q

confidence(compute_χ(events, p_q, p_a))

@test confidence(compute_χ(events, p_q, p_a)) + confidence(compute_χ(events, p_a, p_q)) ≈ 1 atol = 1e-6

##### the only thing I see is a problem in the probabilities themselves?
