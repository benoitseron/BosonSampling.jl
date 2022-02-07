### contains functions to compute probabilities in subsets
### which we used to call partitions
### the functions are in the process of being rewritten for
### multisets

function matrix_phi(k, U, occupancy_vector, n)

	"""matrix P(phi_k) in the notes

	n = number of photons
	occupancy_vector is a Int64 vector of {0,1}^m with one if the mode belongs to the partition, zero otherwise

	remark : opposite sign convention from the notes called notes_singlemode_output
	by default for the r = (1,1, ..., 1) input """

	if k > n
		throw(ArgumentError("k > n"))
	end

	check_at_most_one_particle_per_mode(occupancy_vector)

	m = size(U,1)
	mat = Matrix{eltype(U)}(I, size(U))

	for i in 1:m
		if occupancy_vector[i] == 1
			mat[i,i] = exp(-2im * pi * k / (n+1))
		end
	end

	mat
end

function proba_partition_partial(; U, S, occupancy_vector, input_state, checks = true)

	"""returns a n+1 sized array giving the probability of have [zero, one, ...]
	photons inside the bins given by occupancy_vector for the interferometer
	given by the matrix U (perfectly indistinguishable photons)

	like proba_partition but with partial distinguishability defined through the S
	matrix, with conventions as in Tichy (note that we take the S matrix to be n*n
	while the interferometer has m modes)

	NOTE :
	in the following, we take U to be m*m
	while M is the scattering matrix, as in Tichy, M_ij = U_{d_i}, _{d_j}
	and have n photons, generally in the first n modes
	the distinguishability matrix is defined as in tichy, to be n*n,
	S_{ij} = <phi_{d_i}|phi_{d_j}>
	this is not a problem as it does not depend on the output partition but be aware of it
	"""
	m = size(U,1)
	n = sum(input_state)

	if n == 0
		throw(ArgumentError("number of photons ill defined, zero cannot be computed but trivial"))
	end

	if size(S,1) != n
		throw(ArgumentError("S matrix does not have the size required for this number of photons"))
	end

	if checks
		if !is_unitary(U)
			throw(ArgumentError("U not unitary"))
		elseif !is_a_gram_matrix(S)
			throw(ArgumentError("S not gram"))
		else
			println("Note : checking unitarity of U, that S is gram, slows down significantly probabilities computations")
		end
	end

	function scattering_matrix_amplitudes(U, input_state, occupancy_vector, k)

	    """U tilde in the handwritten notes, corresponds to U^dagger Lambda(eta) U"""

	    scattering_matrix(U' * matrix_phi(k, U, occupancy_vector, sum(input_state)) * U, input_state, input_state)

	end

	proba_fourier = 1/vector_factorial(input_state) .* [ryser(S .* scattering_matrix_amplitudes(U, input_state, occupancy_vector, k)) for k in 0:n]

	if length(proba_fourier) == 0
		throw(ArgumentError("cannot compute the idft of an empty array"))
	end

	p_recovered = [1/(n+1) * sum([exp(2im * pi * j * k / (n+1)) * proba_fourier[k+1] for k in 0:n]) for j in 0:n] # direct inverse fourier transform

	if any([!isapprox(imag(p_recovered[i]),  0., atol = 1e-7, rtol = 1e-5) for i in 1:length(p_recovered)])
		println("WARNING : proba[$j] has a significant imaginary part of $(imag(tmp))")
	elseif !(isapprox(sum(real.(p_recovered)), 1., rtol = 1e-5))
		println("WARNING : probabilites do not sum to one ($(sum(real.(tmp)))), maybe be erreneous")
	end

	real.(p_recovered)

end

function proba_partition_bosonic(;U, occupancy_vector, input_state = ones(Int, size(U,1)), checks = true)

	"""indistinguishable version of proba_partition_partial"""

	n = sum(input_state)

	proba_partition_partial(U = U , S = ones(n,n), occupancy_vector = occupancy_vector, input_state = input_state, checks = checks)

end


function proba_partition_distinguishable(;U, occupancy_vector, input_state = ones(Int, size(U,1)), checks = true)

	"""indistinguishable version of proba_partition_partial"""

	n = sum(input_state)

	proba_partition_partial(U = U , S = Matrix{eltype(U)}(I,n,n), occupancy_vector = occupancy_vector, input_state = input_state, checks = checks)

end
function partition_probability_distribution_distinguishable_rand_walk(part, U)
	"""generates a vector giving the probability to have k photons in
	the partition part for the interferometer U by the random walk method
	discussed in a 22/02/21 email
	part is written like (1,2,4) if it is the output modes 1, 2 and 4"""

	function proba_photon_in_partition(i, part, U)
		"""returns the probability that the photon i ends up in the set of outputs part, for distinguishable particles only"""

		sum([abs(U[i, j])^2 for j in part]) # note the inversion line column
	end


	function walk(probability_vector, walk_number, part, U)

		new_probability_vector = similar(probability_vector)
		n = size(U)[1]
		proba_this_walk = proba_photon_in_partition(walk_number, part, U)

		for i in 1 : n+1
			new_probability_vector[i] = (1-proba_this_walk) * probability_vector[i] + proba_this_walk * probability_vector[i != 1 ? i-1 : n+1]
		end

		new_probability_vector
	end

	n = size(U)[1]
	probability_vector = zeros(n+1)
	probability_vector[1] = 1

	for walk_number in 1 : n
		probability_vector = walk(probability_vector, walk_number, part, U)
	end

	if !(isapprox(sum(probability_vector), 1., rtol = 1e-8))
		println("WARNING : probabilites do not sum to one ($(sum(probability_vector)))")
	end

	probability_vector
end
