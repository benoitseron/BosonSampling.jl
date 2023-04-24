function is_hermitian(M; atol = ATOL)
	isapprox(M, M', atol = atol)
end

function is_positive_semidefinite(M; atol = ATOL)
	eigenvalues = eigvals(M)
	!any(real.(eigenvalues) .< - atol) && !any(abs.(imag.(eigenvalues)) .> atol)
end

function is_a_gram_matrix(M; atol = ATOL)
	is_hermitian(M, atol = atol) && is_positive_semidefinite(M, atol = atol) && all([isapprox(M[i,i], one(eltype(M)), atol = atol) for i in 1:size(M)[1]])
end

function is_unitary(U; atol = ATOL)
    isapprox(U' * U, Matrix(I, size(U)), atol = atol)
end

function is_orthonormal(U; atol = ATOL)
	for i = 1:size(U)[2]
		isapprox(abs(norm(U[:, i])), 1., atol = atol) ? nothing : return false
		for j = 1:size(U)[2]
			if j < i
				isapprox(abs(dot(U[:, i], U[:, j])), 0., atol = atol) ? nothing : return false
			end
		end
	end
end

function is_column_normalized(M; atol = ATOL)

	for i = 1:size(M)[2]
		@test  abs(norm(M[:, i])) ≈ 1. atol = atol
	end
end
