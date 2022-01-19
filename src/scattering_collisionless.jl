# in this file we apply the new tricks of seeing an event probability as
# the probability of finding n photons in the partition defined by the
# output modes, in the case of a collisionless event

include("scattering.jl")

permanent = permanent_ryser


function bunching_probability_brute_force_bosonic(U, input_state, output_state; print_output = false)

    """bosonic bunching probability by direct summation of all possible cases

    bunching_event_proba gives the probability to get the event of [1^n 0^(m-n)]"""

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
