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
@withkw mutable struct QuantumDot <: Source
    efficiency::Real = 1. # probability that a photon is generated if we ask one in this position
end