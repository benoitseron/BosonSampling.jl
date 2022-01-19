# in this file we apply the new tricks of seeing an event probability as
# the probability of finding n photons in the partition defined by the
# output modes, in the case of a collisionless event

include("scattering.jl")

permanent = permanent_ryser



n = 4
m = n^2
input_state = first_modes_array(n,m)
output_state = input_state
partition = output_state

U = rand_haar(m)
S = ones(n,n)

@test is_collisionless(output_state)

H = H_matrix(U,input_state,partition)

permanent(H .* S)

process_probability_partial(U, S, input_state, output_state)
process_probability(U, input_state, output_state)
process_probability_distinguishable(U, input_state, output_state)



bunching_proba, bunching_proba_array, bunching_event_proba = bunching_probability_brute_force_bosonic(U, input_state, output_state, print_output = true)

bunching_event_proba/bunching_proba

function haar_average_bunching_event_proba_versus_bunching_proba(n, m, n_iter)
    input_state = first_modes_array(n,m)
    output_state = input_state
    partition = output_state

    bunching_proba_sum = 0
    bunching_event_proba_sum = 0

    for i in 1:n_iter

        U = rand_haar(m)

        bunching_proba, bunching_proba_array, bunching_event_proba =
        bunching_probability_brute_force_bosonic(U, input_state, output_state, print_output = false)

        bunching_proba_sum+=bunching_proba
        bunching_event_proba_sum+=bunching_event_proba
    end

    bunching_event_proba_sum/bunching_proba_sum
end

haar_average_bunching_event_proba_versus_bunching_proba(3,10,100)
