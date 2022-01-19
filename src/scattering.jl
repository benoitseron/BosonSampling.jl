function iterate_until_collisionless(f)

    """samples f until the result is collisionless

	as it takes a function, example usage :

		random_mode_occupation_collisionless(n::Int,m::Int) = iterate_until_collisionless(() -> random_mode_occupation(n,m))
	"""

    while true

        state = f()
        if is_collisionless(state)
            return state
        end
    end
end

function fill_arrangement(occupation_vector)

    """from a mode occupation list to a mode assignment list (Tichy : s to d(s))

    For instance if vect{s} = (2,0,1) then vect{d} = (1,1,3) """

    arrangement = zeros(eltype(occupation_vector), sum(occupation_vector))

    position = 1
    for i in 1:length(occupation_vector)
        for j in 1:occupation_vector[i]
            arrangement[position] = i
            position += 1
        end
    end

    arrangement
end

fill_arrangement(r::ModeOccupation) = fill_arrangement(r.state)

function random_occupancy(n::Int, m::Int)

	""" returns a vector of size m with n randomly placed ones """
	@warn "swapped argument orders"
	if n > m
		throw(ArgumentError("Implemented at most one photon per mode"))
	else
		occupancy_vector = shuffle(append!(ones(n), zeros(m-n)))
		occupancy_vector = Int.(occupancy_vector)
	end

	return occupancy_vector
end

random_mode_occupation(n::Int,m::Int) = ModeOccupation(random_occupancy(n,m))

function random_mode_occupation_collisionless(n::Int,m::Int)

	n<=m ? iterate_until_collisionless(() -> random_mode_occupation(n,m)) : error("n>m cannot make for collisionless mode occupations")

end

function at_most_one_photon_per_bin(occupancy_vector::Vector{Int})
	all(in([0,1]).(occupancy_vector))
end

check_at_most_one_particle_per_mode(occ) = at_most_one_photon_per_bin(occ) ? error("more than one input per mode") : nothing

function occupancy_vector_to_partition(occupancy_vector)

	"""returns a partition of occupied indexes (mode1,mode2,...) for an occupancy_vector"""
	partition = []

	check_at_most_one_particle_per_mode(occupancy_vector)
	# a partition cannot be defined by taking twice the same mode, it would lead to
	# hidden errors in the results (double counting)

	for (mode, selected) in enumerate(occupancy_vector)
		if selected == 1
			append!(partition, [mode])
		end
	end
	partition

end

function occupancy_vector_to_mode_occupancy(occupancy_vector)

	"""the same as occupancy_vector_to_partition but renamed for clarity when used in the case of, for instance, photon input states"""

	occupancy_vector_to_partition(occupancy_vector)
end


function scattering_matrix(U::Matrix, input_state::Vector{Int}, output_state::Vector{Int})

    """
    U = interferometer matrix, size m*m
    input_state = input occupation number vector (s[i] is the number of photons in the ith input)
    output_state = output occupation number vector (r[i] is the number of photons in the ith output)
    n photons

    follows http://arxiv.org/abs/quant-ph/0406127v1
    """

    m = size(U,1)
    n = sum(input_state)

    if length(input_state) != m || length(output_state) != m
        throw(DimensionMismatch())
    elseif sum(input_state) != sum(output_state)
        error("photon number not the same at input and output")
    end

    index_input = fill_arrangement(input_state)
    index_output = fill_arrangement(output_state)

	U[index_input,index_output]

end

scattering_matrix(interf::Interferometer, r::ModeOccupation, s::ModeOccupation) = scattering_matrix(interf.U,r.state,s.state)
scattering_matrix(interf::Interferometer, i::Input, o::FockDetection) = scattering_matrix(interf.U,i.r,o.s)

function vector_factorial(occupancy_vector)

    """returns the function mu at the denominator of event probabilities"""

    prod(factorial.(occupancy_vector))

end

vector_factorial(r::ModeOccupation) = vector_factorial(r.state)
vector_factorial(i::Input) = vector_factorial(i.r)
vector_factorial(o::FockDetection) = vector_factorial(o.s)

function bosonic_amplitude(U, input_state, output_state, permanent = permanent_ryser)

    """event amplitude"""

    permanent(scattering_matrix(U, input_state, output_state))/sqrt(vector_factorial(input_state) * vector_factorial(output_state))
end

function process_amplitude(U, input_state, output_state, permanent = permanent_ryser)

	bosonic_amplitude(U, input_state, output_state, permanent)
end

function bosonic_probability(U, input_state, output_state)

	"""bosonic process_probability"""

    abs(process_amplitude(U, input_state, output_state))^2

end

function process_probability(U, input_state, output_state)
	#@warn "obsolete function, use bosonic_probability or probability"
	bosonic_probability(U, input_state, output_state)
end

function distinguishable_probability(U, input_state, output_state, permanent = permanent_ryser)

	"""distinguishable (or classical) process_probability"""

	permanent(abs.(scattering_matrix(U, input_state, output_state)).^2)/sqrt(vector_factorial(input_state) * vector_factorial(output_state))

end



function process_probability_distinguishable(U, input_state, output_state, permanent = permanent_ryser)

	#@warn "obsolete function, use distinguishable_probability or probability"

	distinguishable_probability(U, input_state, output_state, permanent)
end


### need to implement partial distinguishability processs probabilities ###

function process_probability_partial(U, S, input_state,output_state)

    """computes the partially distinguishable process probability according to Tichy's tensor permanent https://arxiv.org/abs/1410.7687

    note : naive permanent implementation so really slow"""

    n = size(S,1)

    if sum(input_state) != sum(output_state)
        throw(ArgumentError("particles not conserved"))
    end

    if sum(input_state) != n
        throw(ArgumentError("S matrix doesnt have the same number of photons as the input"))
    end

    M = scattering_matrix(U, input_state, output_state)

    result = zero(eltype(U))

    for sigma in permutations(collect(1:n))
        for rho in permutations(collect(1:n))

            this_diagonal = one(eltype(U))
            for i=1:n
                this_diagonal *= M[sigma[i], i] * conj(M[rho[i],i]) * S[rho[i], sigma[i]]
            end

            result += this_diagonal
        end
    end

    1/(vector_factorial(input_state) * vector_factorial(output_state)) * result

end

process_probability_partial(interf::Interferometer, input_state::Input{PartDist},output_state::OutputMeasurement{FockDetection}) = process_probability_partial(interf.U, input_state.G.S, input_state.r.state,output_state.s.state)

function compute_probability!(ev::Event{TIn,TOut}) where {TIn<:InputType, TOut<:OutputMeasurementType}

	if ev.proba_params.probability != nothing
		@warn "probability was already set in, rewriting"
	end

	if TOut == FockDetection
		if TIn in [Bosonic, Distinguishable, PartDist]

			ev.proba_params.precision = eps()
			ev.proba_params.failure_probability = 0

			if TIn == Bosonic

				ev.proba_params.probability = bosonic_probability(ev.interferometer.U, ev.input_state.r.state, ev.output_measurement.s.state)

			elseif TIn == Distinguishable

				ev.proba_params.probability = distinguishable_probability(ev.interferometer.U, ev.input_state.r.state, ev.output_measurement.s.state)

			else
				ev.proba_params.probability = process_probability_partial(ev.interferometer, ev.input_state, ev.output_measurement)
			end
				ev.proba_params.probability = clean_proba(ev.proba_params.probability)
		end
	else
		error(TOut, " not implemented")
	end
end

function H_matrix(U, input_state, partition_occupancy_vector)

	"""Shshesnovitch's H matrix for a partition defined by partition_occupancy_vector"""

	part = occupancy_vector_to_partition(partition_occupancy_vector)
	input_modes = occupancy_vector_to_mode_occupancy(input_state)

	number_photons = sum(input_state)
	if number_photons == 0
		throw(DomainError("number_photons = 0"))
	end

	if size(U,1) < number_photons
		println("WARNING : more photons than modes (trivially gives the identity as photons are conserved)")
	end

	H = Matrix{ComplexF64}(undef,number_photons,number_photons)

	for i in 1: number_photons
		for j in 1:number_photons
			H[i,j] = sum([U[input_modes[i], l] * conj(U[input_modes[j], l]) for l in part])
		end
	end

	H

end

H_matrix(interf::Interferometer, i::Input, o::OutputMeasurement{FockDetection}) = H_matrix(interf.U, i.r.state, o.s)

### this is an old function that needs to be cleaned

# function partition_probability_distribution_distinguishable(part, U)
# 	"""generates a vector giving the probability to have k photons in
# 	the partition part for the interferometer U by the random walk method
# 	discussed in a 22/02/21 email
# 	part is written like (1,2,4) if it is the output modes 1, 2 and 4"""
#
# 	function proba_photon_in_partition(i, part, U)
# 		"""returns the probability that the photon i ends up in the set of outputs part, for distinguishable particles only"""
#
# 		sum([abs(U[i, j])^2 for j in part]) # note the inversion line column
# 	end
#
#
# 	function walk(probability_vector, walk_number, part, U)
#
# 		new_probability_vector = similar(probability_vector)
# 		n = size(U)[1]
# 		proba_this_walk = proba_photon_in_partition(walk_number, part, U)
#
# 		for i in 1 : n+1
# 			new_probability_vector[i] = (1-proba_this_walk) * probability_vector[i] + proba_this_walk * probability_vector[i != 1 ? i-1 : n+1]
# 		end
#
# 		new_probability_vector
# 	end
#
# 	n = size(U)[1]
# 	probability_vector = zeros(n+1)
# 	probability_vector[1] = 1
#
# 	for walk_number in 1 : n
# 		probability_vector = walk(probability_vector, walk_number, part, U)
# 	end
#
# 	if !(isapprox(sum(probability_vector), 1., rtol = 1e-8))
# 		println("WARNING : probabilites do not sum to one ($(sum(probability_vector)))")
# 	end
#
# 	probability_vector
# end

function is_collisionless(r)
    all(r .<= 1)
end

is_collisionless(r::ModeOccupation) = is_collisionless(r.state)
is_collisionless(i::Input) = is_collisionless(i.r.state)

first_modes_array(n::Int,m::Int) = first_modes(n,m).state
