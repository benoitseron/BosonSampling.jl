function tvd(a,b)
	"""total variation distance"""
	0.5*sum(abs.(a-b))
end

function sqr(a,b)
	"""euclidian distance"""
	sqrt(sum((a-b).^2))
end
