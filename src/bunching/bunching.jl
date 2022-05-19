"""
	H_matrix(U, input_state, partition_occupancy_vector)
	H_matrix(interf::Interferometer, i::Input, o::FockDetection)
	H_matrix(interf::Interferometer, i::Input, subset_modes::ModeOccupation)

H matrix for a partition defined by `partition_occupancy_vector`, see definition
in the article below

**note**: conventions follow the author's [Boson bunching is not
maximized by indistinguishable particles](https://arxiv.org/abs/2203.01306)
which are the ones compatible with Tichy's conventions (Shshnovitch has a
different one for the evolution of the creation operators)
"""
function H_matrix(U, input_state::Vector, partition_occupancy_vector::Vector)

	part = occupancy_vector_to_partition(partition_occupancy_vector)
	input_modes = occupancy_vector_to_mode_occupancy(input_state)

	number_photons = sum(input_state)
	if number_photons == 0
		error("no input photons")
	end

	if size(U,1) < number_photons
		@warn "more photons than modes (trivially gives the identity as photons are conserved)"
	end

	H = Matrix{ComplexF64}(undef,number_photons,number_photons)

	for i in 1: number_photons
		for j in 1:number_photons
			H[i,j] = sum([conj(U[l, input_modes[i]]) * U[l,input_modes[j]] for l in part])
		end
	end

	H

end

H_matrix(interf::Interferometer, i::Input, o::FockDetection) = H_matrix(interf.U, i.r.state, o.s)

H_matrix(interf::Interferometer, i::Input, subset_modes::ModeOccupation) = isa_subset(subset_modes) ? H_matrix(interf.U, i.r.state, subset_modes.state) : error("invalid subset")

"""

	full_bunching_probability(interf::Interferometer, i::Input, subset_modes::Subset)

computes the probability that all n photons end up in the subset of chosen
output modes following [Universality of Generalized Bunching and
Efficient Assessment of Boson Sampling](https://arxiv.org/abs/1509.01561)
"""
function full_bunching_probability(interf::Interferometer, i::Input, subset_modes::Subset)

	return clean_proba(permanent(H_matrix(interf,i,subset_modes) .* transpose(i.G.S)))

end

function full_bunching_probability(interf::Interferometer, i::Input, mo::ModeOccupation)

	return clean_proba(permanent(H_matrix(interf,i,mo) .* transpose(i.G.S)))

end

#
# """
#
# 	bunching_events(input_state::Input, sub::Subset)
#
# generates the output configurations corresponding to a full
# bunching in the subset_modes
# """
# function bunching_events(input_state::Input, sub::Subset)
#
# 	#photon_distribution_in_subset_modes =
# 	all_mode_configurations(input_state, sub, only_photon_number_conserving = false)
#
# 	######### convert to output ModeOccupations
#
# end

"""

	bunching_probability_brute_force_bosonic(U, input_state, output_state; print_output = false)
	bunching_probability_brute_force_bosonic(interf::Interferometer, i::Input, subset_modes::ModeOccupation)

bosonic bunching probability by direct summation of all possible cases

bunching_event_proba gives the probability to get the event of [1^n 0^(m-n)]
"""
function bunching_probability_brute_force_bosonic(U, input_state, output_state; print_output = false)

    n = sum(input_state)
    m = size(U,1)

    bunching_proba = 0
    bunching_proba_array = []
    bunching_event_proba = nothing

    print_output ? println("bunching probabilities : ") : nothing

    for t in reverse.(Iterators.product(fill(0:n,n)...))[:]
        if sum(t) == n # cases where t is physical
            output_state = zeros(Int,m)
            output_state[1:n] .= t[1:n]
            this_proba = process_probability(U, input_state, output_state)
            bunching_proba += this_proba
            push!(bunching_proba_array,[t, this_proba])
            print_output ? println("output = ", t, " p = ", this_proba) : nothing

            if output_state[1:n] == ones(Int, n)
                bunching_event_proba = this_proba
            end
        end
    end

    H = H_matrix(U,input_state,partition)
    @test bunching_proba ≈ real(permanent(H))

    return bunching_proba, bunching_proba_array, bunching_event_proba

end

bunching_probability_brute_force_bosonic(interf::Interferometer, i::Input, subset_modes::ModeOccupation) = bunching_probability_brute_force_bosonic(interf.U, i.r.state, subset_modes.r.state; print_output = false)
