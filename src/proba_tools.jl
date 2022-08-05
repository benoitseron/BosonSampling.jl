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

function remove_nothing(trials)

    new_trials = []

    for trial in trials

        if isa(trial, Number)

            push!(new_trials, trial)

        end

    end

    trials = new_trials

end

"""
    get_power_law_log_log(x_data,y_data)

Gets a power law of type y = exp(c) * x^m from data that looks like a line in a loglog plot.
"""
function get_power_law_log_log(x_data,y_data)

	@argcheck all(x_data .> 0)
	@argcheck all(y_data .> 0)

    x_data = log.(x_data)
    y_data = log.(y_data)

    lr = linregress(x_data,y_data)

    m,c = lr.coeffs

    println("power law: y = $(exp(c)) * x^$m")

    power_law(x) = exp(c)x^(m)

    (power_law, m, c)

end
