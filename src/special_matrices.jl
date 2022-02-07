function fourier_matrix(n; normalized = true)

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

function rand_haar(n)
    """generates a Haar distributed unitary matrix of size n*n"""
    # follows https://case.edu/artsci/math/esmeckes/Meckes_SAMSI_Lecture2.pdf

    qr(randn(ComplexF64, n,n)).Q
end

function matrix_test(n)

    """matrix of 1,2,3... for testing your code"""

    U = Array{Int64}(undef, n, n)

    for i in 1:n
        for j in 1:n
            U[i,j] = (i-1) * n + j
        end
    end

    U
end


function antihermitian_test_matrix(n)
    h = randn(ComplexF64, n,n)

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

function hermitian_test_matrix(n)
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

function permutation_matrix_special_partition_fourier(n)

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

function rand_gram_matrix(n)

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

function rand_gram_matrix_positive(n)

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

	range_rows = collect(1:n)
	range_col = collect(1:n)

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
