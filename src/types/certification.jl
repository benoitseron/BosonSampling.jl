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


mutable struct BayesianPartition <: Certifier
    events::Vector{Event} # input data as events - note that they shouldn't have probabilities associated, just observations
    probabilities::Vector{Real} # array containing the condifence updated at each run for plotting purposes
    confidence::Union{Real, Nothing} # gives the confidence that the null_hypothesis is true
    n_events::Int
    null_hypothesis::Union{HypothesisFunction, InputType}
    alternative_hypothesis::Union{HypothesisFunction, InputType}
    # an input type such as Bosonic can be specified and the correct HypothesisFunction will be created
    part::Union{Partition, Nothing}
    n_subsets::Int

    # function BayesianPartition(events, null_hypothesis::HypothesisFunction, alternative_hypothesis::HypothesisFunction; n_subsets = 2)
    #
    #     if !isa(events, Vector{Event})
    #         events = convert(Vector{Event}, events)
    #     end
    #
    #     for event in events
    #         check_probability_empty(event, resetting_message = false)
    #     end
    #
    #     ev = events[1]
    #     part = equilibrated_partition(ev.input_state.m, n_subsets)
    #     new(events, Vector{Real}(), nothing, length(events), null_hypothesis, alternative_hypothesis, part)
    # end
    #
    # BayesianPartition(events, null_hypothesis::Function, alternative_hypothesis::Function; n_subsets = 2) = Bayesian(events, HypothesisFunction(null_hypothesis), HypothesisFunction(alternative_hypothesis), n_subsets = n_subsets)

    function BayesianPartition(events, null_hypothesis::TIn1, alternative_hypothesis::TIn2, part::Partition) where {TIn1 <: Union{Bosonic, Distinguishable}} where {TIn2 <: Union{Bosonic, Distinguishable}}


        if !isa(events, Vector{Event})
            events = convert(Vector{Event}, events)
        end

        for event in events
            check_probability_empty(event, resetting_message = false)
            ########### need to add a check that always same interferometer, input
        end

        @argcheck TIn1 != TIn2 "no alternative_hypothesis"

        ev = events[1]
        input_modes = ev.input_state.r
        interf = ev.interferometer

        ib = Input{Bosonic}(input_modes)
        id = Input{Distinguishable}(input_modes)

        o = PartitionCountsAll(part)

        evb = Event(ib,o,interf)
        evd = Event(id,o,interf)

        pb = compute_probability!(evb)
        pd = compute_probability!(evd)

        p_partition_B(ev) = p_partition(to_partition_count(ev, part), evb)
        p_partition_D(ev) = p_partition(to_partition_count(ev, part), evd)

        if TIn1 == Bosonic
            p_q = HypothesisFunction(p_partition_B)
            p_a = HypothesisFunction(p_partition_D)
        elseif TIn1 == Distinguishable
            p_a = HypothesisFunction(p_partition_B)
            p_q = HypothesisFunction(p_partition_D)
        end

        new(events, Vector{Real}(), nothing, length(events), p_q, p_a, part, part.n_subset)

    end
end


mutable struct FullBunching <: Certifier
    events::Vector{Event} # input data as events - note that they shouldn't have probabilities associated, just observations
    confidence::Union{Real, Nothing} # gives the confidence that the null_hypothesis is true
    null_hypothesis::Union{HypothesisFunction, InputType}
    alternative_hypothesis::Union{HypothesisFunction, InputType}
    # an input type such as Bosonic can be specified and the correct HypothesisFunction will be created
    subset::Subset
    subset_size::Int

    function FullBunching(events, null_hypothesis::TIn1, alternative_hypothesis::TIn2, subset_size::Int) where {TIn1 <: Union{Bosonic, Distinguishable}} where {TIn2 <: Union{Bosonic, Distinguishable}}

        if !isa(events, Vector{Event})
            events = convert(Vector{Event}, events)
        end

        for event in events
            check_probability_empty(event, resetting_message = false)
            ########### need to add a check that always same interferometer, input
        end

        @argcheck TIn1 != TIn2 "no alternative_hypothesis"

        ev = events[1]
        input_modes = ev.input_state.r

        @argcheck (subset_size > 0 && subset_size < m) "invalid subset"

        @argcheck n*(m-subset_size) > SAFETY_FACTOR_FULL_BUNCHING * m "invalid subset size for high bunching probability"

        subset = Subset(first_modes(subset_size, input_modes.m))

        new(events, nothing, null_hypothesis, alternative_hypothesis,subset, subset_size)

    end
end
