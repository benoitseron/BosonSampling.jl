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

function clean_pdf!(A::Array, atol = ATOL)

	"""checks if an array has all elements as acceptable probabilities within atol
	and converts them to that and summing to one within length(A) * atol
	and renormalizes"""

	A = clean_proba.(A)
	normalization = sum(A)
	if isapprox(normalization, 1, atol = length(A) * atol)
		A = 1/normalization * A
	else
		error("A not normalized")
	end
end
