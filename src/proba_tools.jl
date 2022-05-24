"""
	clean_proba_(probability:Number, atol=ATOL)

Checks whether a (complex) number is close enough to a valid probability with tolerance
`ATOL`. If so, convert it to a positive real number.
"""
function clean_proba(probability::Number, atol=ATOL)

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


"""
	clean_pdf(A::Array, atol=ATOL)

Checks wether an array is an acceptable discrete probability distribution with
tolerance `ATOL`. If so, converts its elements to normalized positive real numbers.
"""
function clean_pdf(A::Array, atol = ATOL)

	"""checks if an array has all elements as acceptable probabilities within atol
	and converts them to that and summing to one within length(A) * atol
	and renormalizes"""

	A .= clean_proba.(A)
	normalization = sum(A)
	if isapprox(normalization, 1, atol = length(A) * atol)
		A .= 1/normalization * A
		return convert(Vector{Real}, A)
	else
		error("A not normalized")
	end
end

"""
	isa_pdf(pdf)

Asserts if `pdf`	is a valid probability distribution.
"""
function isa_pdf(pdf)

	"""asserts if pdf is a probability distribution"""
	clean_pdf(pdf)

end

"""
	tvd(a,b)

Computes the total variation distance between two probability distributions.
"""
function tvd(a,b)
	"""total variation distance"""
	sum(abs.(a-b))
end

"""
	sqr(a,b)

Computes the euclidian distance between two probability distributions.
"""
function sqr(a,b)
	"""euclidian distance"""
	sqrt(sum((a-b).^2))
end
