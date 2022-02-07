include("counter_example_circuit.jl")


println("WARNING this needs to be clarified : U or U'")

input_state = [1 for i in 1:7]
partition_occupancy_vector = [0 for i in 1:7]
partition_occupancy_vector[1] = 1
partition_occupancy_vector[2] = 1

p_pd = proba_partition_partial(U = V, S = S, occupancy_vector = partition_occupancy_vector, input_state = input_state)

p_b = proba_partition_partial(U = V, S = ones(Int64, (7,7)), occupancy_vector = partition_occupancy_vector, input_state = input_state)

scatter([i for i in 0:7], p_pd, label = "partial dist")
scatter!([i for i in 0:7], p_b, label = "bosonic")

### absolute bunching probabilties ###
@test p_pd[8] / p_b[8] > 1.07
@test p_pd[8] ≈  0.007510233224009247
@test p_b[8] ≈ 0.006994170310475809

### distribution of the photons when full bunching ###

n = 7

p_mode_1_pd = zeros(ComplexF64, n+1)
p_mode_1_b = zeros(ComplexF64, n+1)

println("WARNING this needs to be clarified : U or U' etc")

W = conj!(U)

#### this requires the conjugate of U
#### so if the things above are true, see it as just U' transpose, so there would be a transpose too many inside the process_probability_partial compared to the proba_partition_partial
for k in 0:n
    output = [0 for i in 1:7]
    output[1] = k
    output[2] = n-k

    p_mode_1_pd[k+1] = process_probability_partial(W, S, input_state, output)
    p_mode_1_b[k+1] = process_probability_partial(W, ones(Int64, (7,7)), input_state, output)
end

scatter([i for i in 0:7], real.(p_mode_1_pd), label = "partial dist")
scatter!([i for i in 0:7], real.(p_mode_1_b), label = "bosonic")

@test real(sum(p_mode_1_pd) / sum(p_mode_1_b)) > 1.07

sum(p_mode_1_pd)
sum(p_mode_1_b)

####### now as I made changes I clearly need to check whether I need U or U' and where (or U transpose ?)

input_state = [1 for i in 1:7]
output_state = [0 for i in 1:7]
output_state[1] = 7
output_state[2] = 7 - output_state[1]

M = scattering_matrix(U, input_state, output_state)

permanent_ryser(M)
