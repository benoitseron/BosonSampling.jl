abstract type Certifier end

struct Bayesian <: Certifier
    events::Vector{Event} # input data as events - note that they shouldn't have probabilities associated, just observations
    probabilities::Vector{Real} # array containing the condifence updated at each run for plotting purposes
    confidence::Union{Real, Nothing}
    n_events::Int

    function Bayesian(events::Vector{Event})
        for event in events
            check_probability_empty(event, resetting_message = false)
        end
        new(events, Vector{Real}(), nothing, length(events))
    end

end
