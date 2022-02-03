function is_hermitian(M)
	M ≈ M'
end

function is_positive_semidefinite(M, tol = 100*eps())
	!any(eigvals(M) .< - tol)
end

function is_a_gram_matrix(M)

	is_hermitian(M) && is_positive_semidefinite(M) && all([M[i,i] ≈ one(eltype(M)) for i in 1:size(M)[1]])

end

function is_unitary(U)
    U' * U ≈ Matrix(I, size(U))
end

function is_orthonormal(U, atol = 100eps())
	for i = 1:size(U)[2]
		@test  abs(norm(U[:, i])) ≈ 1. atol = atol
		for j = 1:size(U)[2]
			if j < i
				@test  abs(dot(U[:, i], U[:, j])) ≈ 0. atol = atol
			end
		end
	end
end

function is_column_normalized(M, atol = 100eps())

	for i = 1:size(M)[2]
		@test  abs(norm(M[:, i])) ≈ 1. atol = atol
	end
end
