# functions relating to counter examples of the bapat sunder conjecture
# and its physical realization

function cholesky_semi_definite_positive(A)

	"""cholesky decomposition (A = R' * R) for a sdp but not strictly positive definite matrix"""

	# method found in some forum
	X = sqrt(A)

	QR_dec = qr(X)

	R = QR_dec.R

	try
		@test R' * R ≈ A
		return R
	catch err
		throw(Exception(err))
		# println(err)
		# println(R)
		# pretty_table(R)
	end

end

function incorporate_in_a_unitary(X)

	"""incorporates the renormalized matrix X in a double sized unitary through the proof of Lemma 29  of Aaronson Arkipov seminal https://arxiv.org/abs/1011.3245"""

	Y = X / norm(X)

	I_yy = Matrix(I, size(Y)) - Y'*Y

	@test is_hermitian(I_yy)
	@test is_positive_semidefinite(I_yy)

	Z = cholesky_semi_definite_positive(I_yy)

	@test Y'*Y + Z'*Z ≈ Matrix(I, size(Y))

	W = rand(ComplexF64, 2 .* size(Y))

	n = size(Y)[1]
	W[1:n, 1:n] = Y
	W[n+1:2n, 1:n] = Z

	W = modified_gram_schmidt(W)

	@test is_unitary(W)

	W
	
end

function incorporate_in_a_unitary_non_square(X)

	"""same as incorporate_in_a_unitary but for a matrix renormalized X of type (m,n) with m >= n
	generates a minimally sized unitary (ex 9*9 interferometer for the 7*2 M' of the first counter example of drury)"""

	# pad X with zeros

	(m,n) = size(X)

	if m<n
		throw(ArgumentError("m<n not implemented"))
	end

	#padding with zeros
	X_extended = zeros(eltype(X), (m,m))
	X_extended[1:m, 1:n] = X

	#now the same as incorporate_in_a_unitary but adapted to take only the minimally necessary elements
	Y = X_extended / norm(X_extended)

	I_yy = Matrix(I, size(Y)) - Y'*Y

	@test is_hermitian(I_yy)
	@test is_positive_semidefinite(I_yy)

	Z = cholesky_semi_definite_positive(I_yy)

	@test Y'*Y + Z'*Z ≈ Matrix(I, size(Y))

	W = rand(ComplexF64, (m+n, m+n))

	W[1:m, 1:n] = Y[1:m, 1:n] #not taking the padding
	W[m+1:m+n, 1:n] = Z[1:n, 1:n]  #taking only the minimal elements necessary to make the first n columns an orthonormal basis

	W = modified_gram_schmidt(W)

	@test is_unitary(W)

	W

end

function add_columns_to_make_square_unitary(M_dagger)

	"""makes a square matrix U of which the first columns are M_dagger, which has to be unitary by itself"""

	n,r = size(M_dagger)
	U = rand(ComplexF64, (n,n))

	U[1:n,1:r] = M_dagger

	U = modified_gram_schmidt(U)

	is_unitary(U)

	U

end
