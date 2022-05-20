"""
	fourier_matrix(n::Int; normalized=true)

Returns a ``n``-by-``n`` Fourier matrix with optional normalization (`true` by default).
"""
function fourier_matrix(n::Int; normalized = true)

    U = Array{ComplexF64}(undef, n, n)

    for i in 1:n
        for j in 1:n
            U[i,j] = exp((2im * pi / n)* (i-1) * (j-1))
        end
    end

    if normalized
        return 1/sqrt(n) * U
    else
        return U
    end
end

"""
	hadamard_matrix(n::Int; normalized=true)

Returns a ``n``-by-``n`` Hadamard matrix with optional normalization (`true` by default).
"""
function hadamard_matrix(n::Int, normalized = true )

	H = Array{ComplexF64}(undef, n, n)

	for i in 1:n
		for j in 1:n
			bi = map(x->parse(Int,x), split(bitstring(Int8(i-1)),""))
			bj = map(x->parse(Int,x), split(bitstring(UInt8(j-1)),""))
			H[i,j] = (-1)^(sum(bi.&bj))
		end
	end

	if normalized
		return 1/sqrt(n) * H
	else
		return H
	end

end

"""
	sylvester_matrix(p::Int; normalized=true)

Returns a ``2^p``-by-``2^p`` Sylvester matrix with optional normalization
(`true` by default) following [https://arxiv.org/abs/1502.06372](https://arxiv.org/abs/1502.06372).
"""
function sylvester_matrix(p; normalized = true)

    """returns the sylvester matrix of size 2^p"""

    function sylvester_element(i,j,p)

        """returns the element A(i,j) for the sylvester matrix of size 2^p"""
        # following https://arxiv.org/abs/1502.06372
        if p == 0
            return 1
        else
            i_b = digits(i-1, base = 2, pad = p+1)
            j_b = digits(j-1, base = 2, pad = p+1)
            return (-1)^(sum(i_b .* j_b))
        end
    end

    A = Array{ComplexF64}(undef, 2^p, 2^p)
    for i in 1:2^p
        for j in 1:2^p
            A[i,j] = sylvester_element(i,j,p)
        end
    end

    normalized == true ? 1/sqrt(2^p) * A : A

end

"""
	rand_haar(n::Int)

Returns a ``n``-by-``n`` Haar distributed unitary matrix following [https://case.edu/artsci/math/esmeckes/Meckes_SAMSI_Lecture2.pdf](https://case.edu/artsci/math/esmeckes/Meckes_SAMSI_Lecture2.pdf).
"""
function rand_haar(n::Int)

    """generates a Haar distributed unitary matrix of size n*n"""
    # follows https://case.edu/artsci/math/esmeckes/Meckes_SAMSI_Lecture2.pdf

    qr(randn(ComplexF64, n,n)).Q
end

function matrix_test(n::Int)

    """matrix of 1,2,3... for testing your code"""

    U = Array{Int64}(undef, n, n)

    for i in 1:n
        for j in 1:n
            U[i,j] = (i-1) * n + j
        end
    end

    U

end

function antihermitian_test_matrix(n::Int)

    h = randn(ComplexF64, n, n)

    for i = 1:n
        h[i,i] = 1im*imag(h[i,i])
    end

    for i = 1:n
        for j = i+1:n
            h[j,i] = -conj(h[i,j])
        end
    end

    h

end

function hermitian_test_matrix(n::Int)

    h = randn(ComplexF64, n,n)

    for i = 1:n
        h[i,i] = real(h[i,i])
    end

    for i = 1:n
        for j = i+1:n
            h[j,i] = conj(h[i,j])
        end
    end

    h

end

function permutation_matrix_special_partition_fourier(n::Int)

    """P * [1,1,1,1,1,0,0,0,0,0] = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]"""

    P = zeros(Int, n,n)

    for i in 1:2:n
        P[i, Int((i+1)/2)] = 1
    end

    for i in 2:2:n
        P[i, Int(n/2 + i/2)] = 1
    end

    P

end

"""
	rand_gram_matrix_from_orthonormal_basis(n::Int, r::Int)

Returns a ``n``-by-``n`` random Gram matrix that generates orthonormal basis of `r` vectors.
"""
function rand_gram_matrix_from_orthonormal_basis(n,r)

	"""generates a random gram matrix with the property that the
	generating vector basis of r vectors is orthonormal (as is
	strangely the case in Drury's work)"""

	function normalized_random_vector(n)
		v = rand(ComplexF64, n)
		1/norm(v) .* v
	end

	if r >= n
		throw(ArgumentError("need rank < dim to have a non trivial result"))
	end

	generating_vectors = 0.
	#while det(generating_vectors) ≈ 0. #check it's a basis
	    generating_vectors = hcat([normalized_random_vector(n) for i = 1:r]...)
	#### TODO caveat : we don't check for linear independance so there is a chance it doesn't work though highly unlikely
	#end

	generating_vectors = modified_gram_schmidt(generating_vectors)

	is_orthonormal(generating_vectors)

	generating_vectors * generating_vectors'

end


function modified_gram_schmidt(input_vectors)

	"""does the modified gram schmidt (numerically stable) on the columns of the matrix input_vectors"""

	function projector(a, u)
		"""projector of a on u"""
		dot(u,a)/dot(u,u) .* u
	end

	final_vectors = copy(input_vectors)

	for i in 1:size(input_vectors)[2]

		if i == 1
			#normalize first vector
			final_vectors[:,1] /= norm(final_vectors[:, 1])
		else
			for j in 1:i-1
				final_vectors[:, i] -= projector(final_vectors[:, i], final_vectors[:, j])
			end
			final_vectors[:, i] /= norm(final_vectors[:, i])

		end

	end

	is_orthonormal(final_vectors)

	final_vectors

end

"""
	rand_gram_matrix(n::Int)

Returns a full rank ``n``-by-``n`` random Gram matrix.
"""
function rand_gram_matrix(n::Int)

	"""random gram matrix with full rank """

	function normalized_random_vector(n)
		v = rand(ComplexF64, n)
		1/norm(v) .* v
	end

	generating_vectors = 0.
	while det(generating_vectors) ≈ 0. #check it's a basis
	    generating_vectors = hcat([normalized_random_vector(n) for i = 1:n]...)

	end

	generating_vectors' * generating_vectors

end

"""
	rand_gram_matrix_rank(n::Int, r::Int)

Returns a ``n``-by-``n`` random Gram matrix of maximum rank and great likelihood `r`.
"""
function rand_gram_matrix_rank(n, r)

	"""random gram matrix with rank at most r, and with great likelihood r
	"""

	function normalized_random_vector(r)
		v = rand(ComplexF64, r)
		1/norm(v) .* v
	end

	generating_vectors = 0.
	generating_vectors = hcat([normalized_random_vector(r) for i = 1:n]...)

	generating_vectors' * generating_vectors


end

"""
	rand_gram_matrix_real(n::Int)

Returns a real ``n``-by-``n`` random Gram matrix of full rank.
"""
function rand_gram_matrix_real(n)

	"""random gram matrix with full rank and real el"""

	function normalized_random_vector(n)
		v = rand(Float64, n)
		1/norm(v) .* v
	end

	generating_vectors = 0.
	while det(generating_vectors) ≈ 0. #check it's a basis
	    generating_vectors = hcat([normalized_random_vector(n) for i = 1:n]...)

	end

	generating_vectors' * generating_vectors

end

"""
	rand_gram_matrix_positive(n::Int)

Returns a positive elements ``n``-by-``n`` random Gram matrix of full rank.
"""
function rand_gram_matrix_positive(n::Int)

	"""random gram matrix with full rank and positive el"""

	function normalized_random_vector(n)
		v = abs.(rand(Float64, n))
		1/norm(v) .* v
	end

	generating_vectors = 0.
	while det(generating_vectors) ≈ 0. #check it's a basis
	    generating_vectors = hcat([normalized_random_vector(n) for i = 1:n]...)

	end

	generating_vectors' * generating_vectors

end

# functions to clear matrices for better readability of a matrix by a human

function set_small_values_to_zero(x, eps)
	abs(x) < eps ? 0.0 : x
end

function set_small_imag_real_to_zero(x, eps = 1e-7)
	if abs(real(x)) < eps
		x = 0.0 + 1im * imag(x)
	end
	if abs(imag(x)) < eps
		x = real(x) + 0.0im
	end
	x
end

function make_matrix_more_readable(M)
	set_small_imag_real_to_zero.(M)
end

function display_matrix_more_readable(M)
	pretty_table(make_matrix_more_readable(M))
end

####### sub matrices ########

function remove_row_col(A, rows_to_remove, cols_to_remove)

	"""hands back the matrix A without the rows labelled in rows_to_remove, columns in cols_to_remove

	ex :
	rows_to_remove = [1]
	cols_to_remove = [2]

	to have A(1,2) in the notations of Minc"""

	n = size(A)[1]
	m = size(A)[2]

	range_rows = collect(1:n)
	range_col = collect(1:m)

	setdiff!(range_rows, rows_to_remove)
	setdiff!(range_col, cols_to_remove)

	A[range_rows, range_col]

end

### permanent related quantities ###

function sub_permanent_matrix(A)

	n = size(A)[1]

	sub_permanents = similar(A)


	for i in 1 : n
		for j in 1 : n

			sub_permanents[i,j] = ryser(remove_row_col(A, [i], [j]))

		end
	end

	sub_permanents

end

function gram_from_n_r_vectors(M)

	# generate a random matrix of n r-vectors M
	# normalize it to have a gram matrix

	n = size(M)[1]

	for i in  1 : n
		M[:,i] ./= norm(M[:,i])
	end

	M' * M # our rank r n*n gram matrix

end

"""
	gram_matrix_one_param(n::Int, x::Real)

Returns a ``n``-by-``n`` Gram matrix parametrized by the real ``0 ≤ x ≦ 1``.
"""
function gram_matrix_one_param(n::Int, x::Real)

	@argcheck x>=0 && x<= 1

	S = 1.0 * Matrix(I, n, n)
	for i in 1:n
		for j in 1:n
			i!=j ? S[i,j] = x : continue
		end
	end

	S

end

function column_normalize(M)
   for col in 1:size(M,2) #column normalize
	   M[:, col] = M[:, col]./norm(M[:, col])
   end
   M
end

"""
	perturbed_gram_matrix(M, epsilon)

Returns a Gram matrix generated by the columns of `M` which are perturbed by
a Gaussian quantity of variance epsilon once normalized.
"""
function perturbed_gram_matrix(M, epsilon)

    """M defines the set of vectors generating the gram matrix
    each column is a generating vector for the gram matrix
    we perturb them by some random gaussian amount with set variance epsilon once normalized """

    M = column_normalize(M)

    d = Normal(0.0, epsilon)

    perturbation_vector = rand(d,size(M))

    M += perturbation_vector

    M = column_normalize(M)

    M' * M

end

"""
	perturbed_unitary(U, epsilon)

Returns a unitary matrix whose columns are generating vector perturbed by a
random Gaussian quantity with variance epsilon once normalized.
"""
function perturbed_unitary(U, epsilon)

    """U a unitary matrix each column is a generating vector
    we perturb them by some random gaussian amount with set variance epsilon once normalized """

    d = Normal(0.0, epsilon)

    perturbation_vector = rand(d,size(U))

    U += perturbation_vector

    U = modified_gram_schmidt(U)

    U

end
