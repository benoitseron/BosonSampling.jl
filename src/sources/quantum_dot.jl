# simulates the input mixture sent by a Quantum Dot (such as used in our OneLoop experiment)

"""
    abstract type Source 

Holds types simulating realistic experimental sources.
"""
abstract type Source end

"""

    mutable struct QuantumDot <: Source

Type holding a model of a QuantumDot. The point is to simulate its non deterministic photon generation, to know what kind of input is sent through the interferometer.

This is held in through the field `efficiency`, the probability that a photon is generated if we ask one in this position. The probabilities are assumed to be IID.
"""
@with_kw mutable struct QuantumDot <: Source
    efficiency::Real = 1. # probability that a photon is generated if we ask one in this position
    d::Distribution = Bernoulli(efficiency)
end

efficiency = 0.8

source = QuantumDot(efficiency = efficiency)

####### create a Type photon distribution, creating a relation between n_in and n_out?


function sample_imperfect_input(state::Vector{Int}, source::Source)