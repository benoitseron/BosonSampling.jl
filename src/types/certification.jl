abstract type Certifier end

"""
    HypothesisFunction

Stores a function acting as a hypothesis. The function `f` eeds to output the probability associated with an `Event` under that hypothesis. It could be any type of `Event`, such as `FockDetection` or `PartitionOccupancy`.
"""
struct HypothesisFunction
    f::Function
end

"""
    Bayesian(events::Vector{Event}, null_hypothesis::HypothesisFunction, alternative_hypothesis::HypothesisFunction)

Certifies using bayesian testing.
"""
mutable struct Bayesian <: Certifier
    events::Vector{Event} # input data as events - note that they shouldn't have probabilities associated, just observations
    probabilities::Vector{Real} # array containing the condifence updated at each run for plotting purposes
    confidence::Union{Real, Nothing} # gives the confidence that the null_hypothesis is true
    n_events::Int
    null_hypothesis::HypothesisFunction
    alternative_hypothesis::HypothesisFunction

    function Bayesian(events, null_hypothesis::HypothesisFunction, alternative_hypothesis::HypothesisFunction)

        if !isa(events, Vector{Event})
            events = convert(Vector{Event}, events)
        end

        for event in events
            check_probability_empty(event, resetting_message = false)
        end
        new(events, Vector{Real}(), nothing, length(events), null_hypothesis, alternative_hypothesis)
    end

    Bayesian(events, null_hypothesis::Function, alternative_hypothesis::Function) = Bayesian(events, HypothesisFunction(null_hypothesis), HypothesisFunction(alternative_hypothesis))

end


# mutable struct BayesianPartition <: Certifier
#     events::Vector{Event} # input data as events - note that they shouldn't have probabilities associated, just observations
#     probabilities::Vector{Real} # array containing the condifence updated at each run for plotting purposes
#     confidence::Union{Real, Nothing} # gives the confidence that the null_hypothesis is true
#     n_events::Int
#     null_hypothesis::HypothesisFunction
#     alternative_hypothesis::HypothesisFunction
#     part::Partition
#
#     function BayesianPartition(events, null_hypothesis::HypothesisFunction, alternative_hypothesis::HypothesisFunction; n_subsets = 2)
#
#         if !isa(events, Vector{Event})
#             events = convert(Vector{Event}, events)
#         end
#
#         for event in events
#             check_probability_empty(event, resetting_message = false)
#         end
#
#         ev = events[1]
#         part = equilibrated_partition(ev.input_state.m, n_subsets)
#         new(events, Vector{Real}(), nothing, length(events), null_hypothesis, alternative_hypothesis, part)
#     end
#
#     BayesianPartition(events, null_hypothesis::Function, alternative_hypothesis::Function; n_subsets = 2) = Bayesian(events, HypothesisFunction(null_hypothesis), HypothesisFunction(alternative_hypothesis), n_subsets = n_subsets)
#
# end
