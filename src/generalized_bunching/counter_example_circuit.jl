include("counter_example_functions.jl")

omega_5 = exp(2/5 * 1im * pi)
M_dagger_unnormalized = 1/sqrt(2) .* [[sqrt(2) 0]; [0 sqrt(2)]; [1 1]; [1 omega_5]; [1 omega_5^2]; [1 omega_5^3]; [1 omega_5^4]]

H_circuit_unnormalized = M_dagger_unnormalized * M_dagger_unnormalized' # unnormalized meaning that it is a gram matrix though the normalization of H must come from U
S_circuit = copy(conj(H_circuit_unnormalized'))

### starting from the M expression ###

M_dagger = sqrt(2/7) .* M_dagger_unnormalized

M_dagger' * M_dagger

H = M_dagger * M_dagger'
H_circuit = H

S = copy(conj((M_dagger_unnormalized * M_dagger_unnormalized')))

@test violates_bapat_sunder(H,S)

U = add_columns_to_make_square_unitary(M_dagger)

V = U'

# the scattering matrix is manifestly U' and not U

###### beam splitter and special counter example matrices ######

function beam_splitter(transmission_amplitude = sqrt(0.5))

	"""2d beam_splitter matrix, follows the conventions of Leonardo"""

	# |t|^2 is the "transmission probability"

	t = transmission_amplitude
	r = sqrt(1-t^2)

	[[t r]; [-r t]]

end

function beam_splitter_modes(;in_up,in_down,out_up,out_down, transmission_amplitude,n)

	"""beam splitter incorrporated for connecting specific input modes, output modes in a dimension n interferometer"""

	if in_up == in_down || out_up == out_down
		throw(ArgumentError())
	end

	bs = Matrix{ComplexF64}(I,n,n)

	sub_bs = beam_splitter(transmission_amplitude)

	bs[in_up, out_up] = sub_bs[1,1]
	bs[in_down, out_down] = sub_bs[2,2]
	# bs[in_down, out_up] = sub_bs[2,1]
	# bs[in_up, out_down] = sub_bs[1,2]
	bs[in_down, out_up] = sub_bs[1,2] ####### this change is a bit ad hoc and needs to be checked carefully with the conventions
	bs[in_up, out_down] = sub_bs[2,1]


	bs

end

function analytical_counter_example_interferometer(n,r)

	"""analytical counter example interferometer as proposed by Leonardo on 7/9/21 for the case n = 7,r=2

	it is not yet verified whether other cases lead to violations"""

	if r >= n
		throw(ArgumentError())
	end

	complete_interferometer = Matrix{ComplexF64}(I,n,n)

	bottom_dft = Matrix{ComplexF64}(I,n,n)
	bottom_dft[r+1:n, r+1:n] = fourier_matrix(n-r)

	complete_interferometer *= bottom_dft

	for top_mode in 1:r
		complete_interferometer *= beam_splitter_modes(in_up = top_mode, in_down = top_mode+r, out_up = top_mode, out_down = top_mode+r, transmission_amplitude = sqrt(r/n), n = n)
	end

	complete_interferometer

end

### comparison with drury ###

@test violates_bapat_sunder(H_circuit,S_circuit)


U_drury = add_columns_to_make_square_unitary(M_dagger)

# U = U'# the scattering matrix is manifestly U' and not U
println("WARNING this needs to be clarified : U or U'")

analytical_counter_example_interferometer(7,2) .≈ U_drury

U_generated_circuit = analytical_counter_example_interferometer(7,2)

H_generated = H_matrix(U_generated_circuit, [1 for i in 1:7], [i <= 2 ? 1 : 0 for i in 1:7])
H_generated_drury = H_matrix(U_drury, [1 for i in 1:7], [i <= 2 ? 1 : 0 for i in 1:7])


@test H_generated ≈ H_circuit

pretty_table(H_generated)

pretty_table(H_circuit)

violates_bapat_sunder(H_generated, S_circuit)

function violation_ratio(n,k)

	"""outputs the violation_ratio for the full bunching of photons in the generalized
	family of circuits

	drury's example :


		n = 7 # number of photons
		k = 2 # number of auxilliary photons"""


	U = analytical_counter_example_interferometer(n,k)

	H= H_matrix(U, [1 for i in 1:n], [i <= k ? 1 : 0 for i in 1:n])

	pol_states = 1/sqrt(k) * ones(ComplexF64, k,n)

	k > 2 ? (@warn "check this is the wanted polarization pattern") : nothing

	omega = exp(2*pi*1im/(n-k))

	pol_states[1:k,1:k] = Matrix(I,k,k)

	for j in 1:k

		pol_states[j,k+1:end] .= 1/sqrt(k) * [omega^((i-1)*(j-1)) for i in 1:(n-k)]

	end

	pol_states

	S = pol_states' * pol_states

	abs(permanent_ryser(H .* S))/abs(permanent_ryser(H))

end

function lower_bound(n)

	"""lower_bound on the violation ratio, see Nat paper"""

	n/8+1/32 * (n-2)^2/(n-1)

end
