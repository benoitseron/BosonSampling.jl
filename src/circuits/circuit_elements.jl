###### beam splitter and special counter example matrices ######

function beam_splitter(transmission_amplitude = sqrt(0.5))

	"""2d beam_splitter matrix, follows the conventions of Leonardo"""

	# |t|^2 is the "transmission probability"

	t = transmission_amplitude
	r = sqrt(1-t^2)

	bs = [[t -r]; [r t]]

end

function beam_splitter_modes(;in_up, in_down, out_up, out_down, transmission_amplitude, n)

	"""beam splitter incorrporated for connecting specific input modes, output modes in a dimension n interferometer"""

	if in_up == in_down || out_up == out_down
		throw(ArgumentError())
	end

	bs = Matrix{ComplexF64}(I, n, n)

	sub_bs = beam_splitter(transmission_amplitude)

	bs[in_up, out_up] = sub_bs[1,1]
	bs[in_down, out_down] = sub_bs[2,2]
	bs[in_down, out_up] = sub_bs[2,1]
	bs[in_up, out_down] = sub_bs[1,2]
	#
	# # before convention checks:
	# bs[in_down, out_up] = sub_bs[1,2] ####### this change is a bit ad hoc and needs to be checked carefully with the conventions
	# bs[in_up, out_down] = sub_bs[2,1]
	#

	bs

end

function rotation_matrix(angle::Float64)
	[cos(angle) -sin(angle);
	 sin(angle) coss(angle)]
 end

 function rotation_matrix_modes(;in_up, in_dow, out_up, out_down, angle, n)

	 """beam splitter incorrporated for connecting specific input modes, output modes in a dimension n interferometer"""

 	if in_up == in_down || out_up == out_down
 		throw(ArgumentError())
 	end

 	r = Matrix{ComplexF64}(I, n, n)

 	sub_r = rotation_matrix(angle)

 	r[in_up, out_up] = sub_r[1,1]
 	r[in_down, out_down] = sub_r[2,2]
 	# bs[in_down, out_up] = sub_bs[2,1]
 	# bs[in_up, out_down] = sub_bs[1,2]
 	r[in_down, out_up] = sub_r[1,2] ####### this change is a bit ad hoc and needs to be checked carefully with the conventions
 	r[in_up, out_down] = sub_r[2,1]

 	r

 end

 function phase_shift(shifted_modes::Array, param_::Array)

	U = zeros(length(shifted_modes), length(shifted_modes))

 	for i in length(shifted_modes)
		for j in length(shifted_modes)
			i == j && i in shifted_modes ? U[i,j] = exp(1im * param_[i]) : continue
		end
	end

	  U

end

# function build_circuit(mode_occ::ModeOccupation, circuit_elements::Vector{Vector{T}}, layers::Vector) where T <: Interferometer
#
# 	m = mode_occ.m
# 	id = Matrix{ComplexF64}(LinearAlgebra.I, m, m)
# 	res = id
#
# 	for i in 1:length(circuit_elements)
# 		u = id
# 		for j in 1:length(circuit_elements[i])
# 			U = circuit_elements[i][j].U
# 			u[layers[i][j], layers[i][j]] = U
# 		end
#
# 		res *= u
# 	end
#
# 	res
# end
#
# mode_occ = first_modes(3,5)
#
# m = mode_occ.m
#
# b1 = BeamSplitter(1/sqrt(6))
#
# is_unitary(b1.U)
#
# U1 = RandHaar(1)
#
# is_unitary(U1.U)
#
# b2 = BeamSplitter(1/sqrt(2))
#
# is_unitary(b2.U)
#
# U2 = RandHaar(1)
#
# is_unitary(U2.U)
#
# circuit_elem = [[b1,U1],[b2,U2]]
#
# v_ = [[[1,2],[3]], [[1,3],[2]]]
#
# res = build_circuit(mode_occ, circuit_elem, v_)
#
# is_unitary(res)
#
# res * conj(transpose(res))
#
