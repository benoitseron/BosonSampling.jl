###### beam splitter and special counter example matrices ######

function beam_splitter(transmission_amplitude = sqrt(0.5))

	"""2d beam_splitter matrix, follows the conventions of Leonardo"""

	# |t|^2 is the "transmission probability"

	t = transmission_amplitude
	r = sqrt(1-t^2)

	bs = [[t -r]; [r t]]

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
	bs[in_down, out_up] = sub_bs[2,1]
	bs[in_up, out_down] = sub_bs[1,2]
	#
	# # before convention checks:
	# bs[in_down, out_up] = sub_bs[1,2] ####### this change is a bit ad hoc and needs to be checked carefully with the conventions
	# bs[in_up, out_down] = sub_bs[2,1]
	#

	bs

end
