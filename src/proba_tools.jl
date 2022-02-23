function clean_proba(probability, atol=ATOL)

    """checks if a (complex) number is a close enough probability
    converts it to a positive real number if complex or just positive if real"""

    not_a_proba() = error(probability , " is not a probability")

    if real(probability) >= -ATOL && real(probability) <= 1 + ATOL
        if isa(probability,Complex)
            if abs(imag(probability)) <= ATOL
                return abs(probability)
            else
                not_a_proba()
            end
        elseif isa(probability,Real)
            return abs(probability)
        else
            error(typeof(probability), " probability type undefined")
        end
    else
        not_a_proba()
    end
end

function tvd(a,b)
	"""total variation distance"""
	0.5*sum(abs.(a-b))
end

function sqr(a,b)
	"""euclidian distance"""
	sqrt(sum((a-b).^2))
end
