function violates_bapat_sunder(A,B, tol = ATOL)

	diagonal_el_product = prod([B[i,i] for i in 1:size(B)[1]])

	# @test (imag(diagonal_el_product) ≈ 0.)


	if !(imag(diagonal_el_product) ≈ 0.)
		throw(Exception("product of diagonal elements of B is not real"))
	elseif !(is_positive_semidefinite(A) && is_positive_semidefinite(B))
		throw(Exception("A or B not semidefinitepositive"))
	else
		return abs(permanent_ryser(A .* B)) / (abs(permanent_ryser(A)) * real(prod([B[i,i] for i in 1:size(B)[1]]))) > 1+tol
	end
	return nothing
end
