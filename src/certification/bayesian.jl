### first, a simple bayesian estimator ###


confidence(χ) = χ/(1+χ)

function update_confidence(event, p_q, p_a, χ)

    χ *= p_q(event)/p_a(event)
    χ

end

function compute_confidence(events,p_q, p_a)

    """a bayesian confidence estimator:

    returns the probability that the null hypothesis Q is right compared to
    the alternative hypothesis A

    the functions take in events
    p_q is a function that takes in an event and gives its probability in the null hypothesis
    p_a is the same for the alternative hypothesis
    χ is the ratio between null hypothesis probabilities and alternative ones

    """

    χ = 1
    for event in events
        χ = update_confidence(event, p_q, p_a, χ)
    end
    confidence(χ)
end
