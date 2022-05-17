"""
	violates_bapat_sunder(A,B, tol = ATOL)

Checks if matrices `A` and `B` violate the Bapat-Sunder conjecture, see
[Boson bunching is not maximized by indistinguishable particles](https://arxiv.org/abs/2203.01306)
"""
function violates_bapat_sunder(A,B, tol = ATOL)

	diagonal_el_product = prod([B[i,i] for i in 1:size(B)[1]])

	# @test (imag(diagonal_el_product) ≈ 0.)

	if !(imag(diagonal_el_product) ≈ 0.)
		throw(Exception("product of diagonal elements of B is not real"))
	elseif !(is_positive_semidefinite(A) && is_positive_semidefinite(B))
		throw(Exception("A or B not semidefinitepositive"))
	else
		return abs(ryser(A .* B)) / (abs(ryser(A)) * real(prod([B[i,i] for i in 1:size(B)[1]]))) > 1+tol
	end

	return nothing

end
