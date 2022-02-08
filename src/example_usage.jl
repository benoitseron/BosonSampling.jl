### scattering ###

m = 5
n = 3

interf = RandHaar(m)
i = Input{PartDist}(first_modes(n,m), GramMatrix{PartDist}(n, rand_gram_matrix(n)))
o = OutputMeasurement{FockDetection}(first_modes(n,m))
ev = Event(i,o,interf)

compute_probability!(ev)

### bunching ###

m = 5
n = 3

interf = RandHaar(m)
ib = Input{Bosonic}(first_modes(n,m))
ipd = Input{PartDist}(first_modes(n,m), GramMatrix{PartDist}(n, rand_gram_matrix(n)))
subset_modes = first_modes(n,m)

pb = full_bunching_probability(interf, ib, subset_modes)
ppd = full_bunching_probability(interf, ipd, subset_modes)

@test pb/ppd > 1. # this doesn't HAVE TO pass but will pass in nearly all
# cases

### classical sampling ###

classical_sampler(U = rand_haar(16), m = 16, n = 3)

### MIS sampling ###

n = 8
m = n^2
starting_state = zeros(Int, m)
input_state = first_modes_array(n,m)

U = copy(rand_haar(m))

# generate a collisionless state as a starting point
starting_state = iterate_until_collisionless(() -> random_occupancy(n,m))

known_pdf(state) = process_probability_distinguishable(U, input_state, state)
target_pdf(state) = process_probability(U, input_state, state)
known_sampler = () -> iterate_until_collisionless(() -> classical_sampler(U = U, m = m, n = n)) # gives a classical sampler


samples = metropolis_sampler(;target_pdf = target_pdf, known_pdf = known_pdf , known_sampler = known_sampler , starting_state = starting_state, n_iter = 100)