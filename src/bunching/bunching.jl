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

H_matrix(interf::Interferometer, i::Input, subset_modes::ModeOccupation) = isa_subset(subset_modes) ? H_matrix(interf.U, i.r.state, subset_modes.state) : error("invalid subset")

function full_bunching_probability(interf::Interferometer, i::Input, subset_modes::ModeOccupation)

	"""computes the probability that all n photons end up in the subset of chosen
	output modes following https://arxiv.org/abs/1509.01561"""

	return clean_proba(permanent(H_matrix(interf,i,subset_modes) .* i.G.S))

end
