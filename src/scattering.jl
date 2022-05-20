"""
	iterate_until_collisionless(f)

Sample `f` until the result is collisionless.
"""
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

"""
	fill_arrangement(occupation_vector)
	fill_arrangement(r::ModeOccupation)
	fill_arrangement(input::Input)

Convert a mode occupation list to a mode assignement.
"""
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

fill_arrangement(inp::Input) = fill_arrangement(inp.r)

"""
	random_occupancy(n::Int, m::Int)

Return a vector of size `m`	with `n` randomly placed ones.
"""
function random_occupancy(n::Int, m::Int)

	""" returns a vector of size m with n randomly placed ones """
	if n > m
		throw(ArgumentError("Implemented at most one photon per mode"))
	else
		occupancy_vector = shuffle(append!(ones(n), zeros(m-n)))
		occupancy_vector = Int.(occupancy_vector)
	end

	return occupancy_vector
end

"""
	random_mode_occupation(n::Int, m::Int)

Create a [`ModeOccupation`](@ref) from a mode occupation list of `n` ramdomly placed ones
among `m` sites.
"""
random_mode_occupation(n::Int, m::Int) = ModeOccupation(random_occupancy(n,m))


"""
	random_mode_occupation_collisionless(n::Int, m::Int)

Create a [`ModeOccupation`](@ref) from a random mode occupation that is likely collisionless.
"""
function random_mode_occupation_collisionless(n::Int, m::Int)

	n<=m ? iterate_until_collisionless(() -> random_mode_occupation(n,m)) : error("n>m cannot make for collisionless mode occupations")

end

"""
	at_most_one_photon_per_bin(occupancy_vector::Vector{Int})
	check_at_most_one_particle_per_mode(occ)

Check wether `occupancy_vector` contains more than one photon per site.
"""
function at_most_one_photon_per_bin(occupancy_vector::Vector{Int})
	all(in([0,1]).(occupancy_vector))
end

check_at_most_one_particle_per_mode(occ) = at_most_one_photon_per_bin(occ) ? nothing : error("more than one input per mode")

"""
	occupancy_vector_to_partition(occupancy_vector)
	occupancy_vector_to_mode_occupancy(occupancy_vector)

Return a partition of occupied modes from an `occupancy_vector`.
"""
function occupancy_vector_to_partition(occupancy_vector)

	"""returns a partition of occupied indexes (mode1,mode2,...) for an occupancy_vector"""
	partition = []

	check_at_most_one_particle_per_mode(occupancy_vector)
	# a partition cannot be defined by taking twice the same mode, it would lead to
	# hidden errors in the results (double counting)

	for (mode, selected) in enumerate(occupancy_vector)
		if selected == 1
			append!(partition, mode)
		end
	end
	partition

end

function occupancy_vector_to_mode_occupancy(occupancy_vector)

	"""the same as occupancy_vector_to_partition but renamed for clarity when used in the case of, for instance, photon input states"""

	occupancy_vector_to_partition(occupancy_vector)
end

"""
	scattering_matrix(U::Matrix, input_state::Vector{Int}, output_state::Vector{Int})
	scattering_matrix(U::Interferometer, r::ModeOccupation, s::ModeOccupation)
	scattering_matrix(U::Interferometer, i::Input, o::FockDetection)

Return the submatrix of `U` whose rows and columns are respectively defined by
`input_state` and `output_state`.

!!! note "Reference"
	[http://arxiv.org/abs/quant-ph/0406127v1](http://arxiv.org/abs/quant-ph/0406127v1)
"""
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

	try
		U[index_input,index_output]
	catch err
		if isa(MethodError, err)
			copy(U)[index_input,index_output]
		else
			println(err)
		end
	end

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

"""
	bosonic_amplitude(U, input_state, output_state, permanent=ryser)
	process_amplitude(U, input_state, output_state, permanent=ryser)

Compute the probability amplitude to go from `input_state` to `output_state`
through the interferomter `U` in the [`Bosonic`](@ref) case.
"""
function bosonic_amplitude(U, input_state, output_state, permanent = ryser)

    """event amplitude"""

    permanent(scattering_matrix(U, input_state, output_state))/sqrt(vector_factorial(input_state) * vector_factorial(output_state))
end

function process_amplitude(U, input_state, output_state, permanent = ryser)

	bosonic_amplitude(U, input_state, output_state, permanent)
end

"""
	bosonic_probability(U, input_state, output_state)
	process_probability(U, input_state, output_state)

Compute the probability to go from `input_state` to `output_state`
through the interferometer `U` in the [`Bosonic`](@ref) case.
"""
function bosonic_probability(U, input_state, output_state)

	"""bosonic process_probability"""

    abs(process_amplitude(U, input_state, output_state))^2

end

function process_probability(U, input_state, output_state)
	#@warn "obsolete function, use bosonic_probability or probability"
	bosonic_probability(U, input_state, output_state)
end

"""
	distinguishable_probability(U, input_state, output_state, permanent=ryser)
	process_probability_distinguishable(U, input_state, output_state, permanent=ryser)

Compute the probability to go from `input_state` to `output_state` through
the interferomter `U` in the [`Distinguishable`](@ref) case.
"""
function distinguishable_probability(U, input_state, output_state, permanent = ryser)

	"""distinguishable (or classical) process_probability"""

	permanent(abs.(scattering_matrix(U, input_state, output_state)).^2)/sqrt(vector_factorial(input_state) * vector_factorial(output_state))

end

function process_probability_distinguishable(U, input_state, output_state, permanent = ryser)

	#@warn "obsolete function, use distinguishable_probability or probability"

	distinguishable_probability(U, input_state, output_state, permanent)
end


### need to implement partial distinguishability processs probabilities ###

"""
	process_probability_partial(U, S, input_state, output_state)
	process_probability_partial(interf::Interferometer, input_state::Input{TIn} where {TIn<:PartDist},output_state::FockDetection)

Compute the probability to go from `input_state` to `output_state` through the
interferometer `U` in the [`PartDist`](@ref) case where partial distinguishable is described
by the [`GramMatrix`](@ref) `S`.

!!! note "Reference"
    [https://arxiv.org/abs/1410.7687](https://arxiv.org/abs/1410.7687)
"""
function process_probability_partial(U, S, input_state,output_state)

    """computes the partially distinguishable process probability according to Tichy's tensor permanent https://arxiv.org/abs/1410.7687"""

    n = size(S,1)

    @argcheck sum(input_state) == sum(output_state) "particles not conserved"

    @argcheck sum(input_state) == n "S matrix doesnt have the same number of photons as the input"

	M = scattering_matrix(U, input_state, output_state)
	W = Array{eltype(U)}(undef, (n,n,n))

	for ss in 1:n
	    for rr in 1:n
	        for j in 1:n
	            W[ss,rr,j] = M[ss, j] * conj(M[rr, j]) * S[rr, ss]
	        end
	    end
	end

    1/(vector_factorial(input_state) * vector_factorial(output_state)) * ryser_tensor(W)

end

process_probability_partial(interf::Interferometer, input_state::Input{TIn} where {TIn<:PartDist},output_state::FockDetection) = process_probability_partial(interf.U, input_state.G.S, input_state.r.state,output_state.s.state)

function compute_probability!(ev::Event{TIn,TOut}) where {TIn<:InputType, TOut<:FockDetection}

	check_probability_empty(ev)

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

"""
	output_mode_occupation(n::Int, m::Int)

Return all possible configurations of `n` photons among `m` modes.
"""
function output_mode_occupation(number_photons, number_modes)

	nlist = collect(1:number_modes)

	for i = 1:number_photons-1
		sublist = []
		for j = 1:length(nlist)
			for k = 1:number_modes
				push!(sublist, vcat(nlist[j], k))
			end
		end
		nlist = sublist
	end

	nlist

end

"""
	check_suppression_law(event)

Check if the event is suppressed according to the [rule](https://arxiv.org/pdf/1002.5038.pdf).
"""
function check_suppression_law(event)

	if mod(sum(event), length(event)) != 0
		return true
	else
		return false
	end

end


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
