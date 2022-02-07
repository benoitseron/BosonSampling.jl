using Plots
using PyPlot
#pyplot()

#Plots.PyPlotBackend()

###### preludes ######

include("counter_example_functions.jl")
include("counter_example_larger_violation_circuits.jl")


omega_5 = exp(2/5 * 1im * pi)
M_dagger_circuit = 1/sqrt(2) .* [[sqrt(2) 0]; [0 sqrt(2)]; [1 1]; [1 omega_5]; [1 omega_5^2]; [1 omega_5^3]; [1 omega_5^4]]

H_circuit_unnormalized = M_dagger_circuit * M_dagger_circuit' # unnormalized meaning that it is a gram matrix though the normalization of H must come from U
S_circuit = copy(conj(H_circuit_unnormalized'))

### starting from the M expression ###

M_dagger = sqrt(2/7) .* M_dagger_circuit
M_dagger' * M_dagger
H = M_dagger * M_dagger'
S = copy(conj((M_dagger_circuit * M_dagger_circuit')))

@test violates_bapat_sunder(H,S)

###### bunching probability ######

U = add_columns_to_make_square_unitary(M_dagger)

U = U'# the scattering matrix is manifestly U' and not U
println("WARNING this needs to be clarified : U or U'")

input_state = [1 for i in 1:7]
partition_occupancy_vector = [0 for i in 1:7]
partition_occupancy_vector[1] = 1
partition_occupancy_vector[2] = 1

p_pd = proba_partition_partial(U = U, S = S, occupancy_vector = partition_occupancy_vector, input_state = input_state)

p_b = proba_partition_partial(U = U, S = ones(Int64, (7,7)), occupancy_vector = partition_occupancy_vector, input_state = input_state)

plt = scatter([i for i in 0:7], p_pd, label = "partially dist", dpi = 600);
scatter!([i for i in 0:7], p_b, label = "bosonic");
plot!(legend = false)
savefig(plt, "images/nature_plots/bunching_probability.png")

@test p_pd[8] / p_b[8] > 1.07


###### mode distribution ######

n = 7

p_mode_1_pd = zeros(ComplexF64, n+1)
p_mode_1_b = zeros(ComplexF64, n+1)

U = copy(conj.(U'))
println("WARNING this needs to be clarified : U or U' etc")

for k in 0:n
    output = [0 for i in 1:7]
    output[1] = k
    output[2] = n-k

    p_mode_1_pd[k+1] = process_probability_partial(U, S, input_state, output)
    p_mode_1_b[k+1] = process_probability_partial(U, ones(Int64, (7,7)), input_state, output)
end

plt = scatter([i for i in 0:7], real.(p_mode_1_pd), label = "partial dist", dpi = 600);
scatter!([i for i in 0:7], real.(p_mode_1_b), label = "bosonic");
plot!(legend = false)
savefig(plt, "images/nature_plots/mode_distribution.png")


@test real(sum(p_mode_1_pd) / sum(p_mode_1_b)) > 1.07

sum(p_mode_1_pd)
sum(p_mode_1_b)


###### mode distribution output of the DFT ######

q_min = 5
q_max = 15
q_array = q_min:q_max

function star_state(q)

	"""star state defined in the nature paper"""

	states = Matrix{ComplexF64}(undef, 2,q)
	states[1,:] .= 1/sqrt(2)
	for j in 1:q
		states[2,j] = 1/sqrt(2) * exp(2pi*1im*j/q)
	end


	states

end

function gram_star_state(q)

	"""corresponding gram matrix"""

	star_state(q)' * star_state(q)
end

function top_two_mode_distribution(U, S)

	"""probability array to get k,n-k,0,...,0 photons as the output mode

	outputs that of the bosonic, distinguishable, and the given S matrix"""

	n = size(U,1)

	p_mode_1_b = zeros(ComplexF64,n+1)
	p_mode_1_pd = copy(p_mode_1_b)
	p_mode_1_d = copy(p_mode_1_b)


	input_state = [1 for i in 1:n]
	partition_occupancy_vector = [0 for i in 1:n]
	partition_occupancy_vector[1] = 1
	partition_occupancy_vector[2] = 1

	for k in 0:n
	    output = [0 for i in 1:n]
	    output[1] = k
	    output[2] = n-k

	    p_mode_1_pd[k+1] = process_probability_partial(U, S, input_state, output)
	    p_mode_1_b[k+1] = process_probability_partial(U, ones(Int64, (n,n)), input_state, output)
		p_mode_1_d[k+1] = distinguishable_probability(U,  input_state, output)
	end

	real.([p_mode_1_b, p_mode_1_pd, p_mode_1_d])
end


function top_two_mode_distribution_renormalized(U, S)

	"""each pdf is normalized to one"""
	return_array = top_two_mode_distribution(U, S)
	for (i,dist) in enumerate(return_array)
		return_array[i] = return_array[i]./sum(return_array[i])
	end
	return_array
end

function plot_top_two_mode_distribution_star_state_normalized(n)

	dists = top_two_mode_distribution_renormalized(fourier_matrix(n), gram_star_state(n))

	plt = plot()

	scatter!([0:n], dists[1], label = "bosonic");
	scatter!([0:n], dists[2], label = "distinguishable");
	scatter!([0:n], dists[3], label = "star state");
	#title!("probability to have k photons in the first mode upon vacuum in the last $(n-2) of a DFT");

	plt
end

plt_array = []

@showprogress for q in q_array
	this_plt = plot_top_two_mode_distribution_star_state_normalized(q)
	push!(plt_array, this_plt)
	savefig(this_plt, "images/nature_plots/probability_to_have_k_photons_in_the_first_mode_upon_vacuum_in_the_last_$(q)_of_a_DFT.png")
	display(this_plt)
end

###### bunching probability compared to lower bound ######

n_bunching = [i for i in 3:12]
bunching_ratio = [0.75, 0.75, 0.8125, 0.933333, 1.07378, 1.22366, 1.37848, 1.53623, 1.69585, 1.85678, 2.01865, 2.18124, 2.34439][1:length(n_bunching)]
lower_bound = [0.5, 0.583333, 0.6875, 0.8, 0.916667, 1.03571, 1.15625, 1.27778, 1.4, 1.52273, 1.64583, 1.76923, 1.89286][1:length(n_bunching)]

plt = scatter(n_bunching, bunching_ratio, dpi = 600);
scatter!(n_bunching, lower_bound);
plot!(legend = false)
hline!([1], linestyle=:dash)
savefig(plt, "images/nature_plots/bunching_ratio.png")

###### interpolating model ######

function interpolation_S_matrix(x; S_violation = S)

	"""a convex combination interpolation between the S matrix of the violation example (x = 0)
	and that of the indistinguishable case"""

	S_indist = ones(size(S_violation))

	(1-x) * S_violation + x * S_indist

end


function proba_n_S_minus_bosonic(Sx; H = H, S = S_circuit)

	"""this gives the probability to get n particles in the partition for
	a S matrix minus that of the indistinguishable case"""

	S_indist = ones(size(S))

	abs(permanent_ryser(H .* Sx)) - abs(permanent_ryser(H .* S_indist))

end

function relative_proba_n_S_minus_bosonic(Sx; H = H, S = S_circuit)

	"""this gives the difference in probability to get n particles in the partition for
	a S matrix minus that of the indistinguishable case relative to the probability of the violation distinguishability matrix"""

	S_indist = ones(size(S))

	(abs(permanent_ryser(H .* Sx)) - abs(permanent_ryser(H .* S_indist)))/abs(permanent_ryser(H .* S))

end
#
# # plot_relative_proba_n_S_minus_bosonic - non minimal example
#
# x_array = range(0,1, length = 101)
# proba_array = relative_proba_n_S_minus_bosonic.(interpolation_S_matrix.(x_array))
#
# plt_interp = plot(x_array, proba_array, label = false)
# xlabel!("x")
# ylabel!("relative violation")

# plot_relative_proba_n_S_minus_bosonic - minimal counter example

x_array = range(0,1, length = 101)
proba_array = relative_proba_n_S_minus_bosonic.(interpolation_S_matrix.(x_array), H = H)

plt_interp = plot(x_array, proba_array, label = false, dpi = 600)
# xlabel!("x")
# ylabel!("relative violation")

savefig(plt_interp, "images/nature_plots/interpolation.png")


##### the case in which p_n^B = p_n^S(x) #####

function proba_n_S_minus_bosonic_interpolation(x)

	proba_n_S_minus_bosonic(interpolation_S_matrix(x))

end

find_zero(proba_n_S_minus_bosonic_interpolation, 0.09)



##### ternary plot #####
#
# see ternary_plot.jl for the data
# as well as the python files for the plotting


###### resiliency to perturbations ######

using Distributions
using Random

include("counter_example_circuit.jl")

H_circuit = H

omega_5 = exp(2/5 * 1im * pi)
M_dagger_circuit = 1/sqrt(2) .* [[sqrt(2) 0]; [0 sqrt(2)]; [1 1]; [1 omega_5]; [1 omega_5^2]; [1 omega_5^3]; [1 omega_5^4]]

M_circuit = conj(copy(M_dagger_circuit')) # the presence of the conjugate is not perfectly clear

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

perturbed_gram_matrix(M_circuit, 0.001)
S_circuit

### plot perturbations ###

function relative_violation(H,S)

	S_bosonic = ones(size(S))

	abs(permanent_ryser(H .* S))/abs(permanent_ryser(H .* S_bosonic))

end

@test relative_violation(H_circuit,S_circuit) > 1.07

@test relative_violation(H_circuit,perturbed_gram_matrix(M_circuit,0)) > 1.07

@test perturbed_gram_matrix(M_circuit,0) ≈ S_circuit

n_trials = 10000
epsilon_array = range(0,0.15, length = 11)

results = zeros(length(epsilon_array), n_trials)

for (i, epsilon) in enumerate(epsilon_array)
	for trial in 1:n_trials
		results[i, trial] = relative_violation(H_circuit,perturbed_gram_matrix(M_circuit,epsilon))
	end
end

#scatter(epsilon_array, results, legend = false)

mean_array = similar(epsilon_array)
var_array = similar(epsilon_array)

for i in 1:length(mean_array)
	mean_array[i] = mean(results[i,:])
	var_array[i] = var(results[i,:])
end

plt_rand = scatter(epsilon_array, mean_array, yerr = sqrt.(var_array), label = "internal states");
hline!([1], linestyle=:dash, label = false);
# xlabel!("epsilon")
# ylabel!("relative violation")


savefig(plt_rand, "images/nature_plots/gram_perturbations.png")

epsilon_array[end-1]



###### random perturbation of matrix elements ######

function perturbed_unitary(U, epsilon)

    """U a unitary matrix each column is a generating vector
    we perturb them by some random gaussian amount with set variance epsilon once normalized """

    d = Normal(0.0, epsilon)

    perturbation_vector = rand(d,size(U))

    U += perturbation_vector

    U = modified_gram_schmidt(U)

    U
end


U_circuit = analytical_counter_example_interferometer(7,2)

U_to_H_circuit(U) = H_matrix(U, [1 for i in 1:7], [i <= 2 ? 1 : 0 for i in 1:7])

@test relative_violation(U_to_H_circuit(U_circuit),S_circuit) > 1.07

@test relative_violation(U_to_H_circuit(U_circuit),perturbed_gram_matrix(M_circuit,0)) > 1.07

@test perturbed_gram_matrix(M_circuit,0) ≈ S_circuit

n_trials = 10000
epsilon_array = range(0,0.09, length = 16)

results = zeros(length(epsilon_array), n_trials)

for (i, epsilon) in enumerate(epsilon_array)
	for trial in 1:n_trials
		results[i, trial] = relative_violation(U_to_H_circuit(perturbed_unitary(U_circuit, epsilon)), S_circuit)
	end
end

#scatter(epsilon_array, results, legend = false)

mean_array_unit = similar(epsilon_array)
var_array_unit = similar(epsilon_array)

for i in 1:length(mean_array)
	mean_array_unit[i] = mean(results[i,:])
	var_array_unit[i] = var(results[i,:])
end

scatter!(plt_rand, epsilon_array, mean_array_unit, yerr = sqrt.(var_array_unit), label = "matrix elements");
ylims!((0.90,1.1));
xlims!((0,0.16))

# xlabel!("epsilon")
# ylabel!("relative violation")

savefig(plt_rand, "images/nature_plots/perturbations.png")
