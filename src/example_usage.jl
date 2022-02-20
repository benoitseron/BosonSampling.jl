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

### subsets ###


s1 = Subset([1,1,0,0,0])
s2 = Subset([0,0,1,1,0])
s3 = Subset([1,0,1,0,0])

"subsets are not allowed to overlap"

check_subset_overlap([s1,s2,s3])

### HOM tests: one mode ###

m = 2
n = 2

set1 = [1,0]
physical_interferometer = Fourier(m)
part = Partition([Subset(set1)])

(physical_indexes,  pdf) = compute_probabilities_partition(physical_interferometer, part, n)

photon_number_conserving_events(physical_indexes,n)

check_photon_conservation(physical_indexes, pdf, n)

print_pdfs(physical_indexes,  pdf, n)

### HOM tests: mode1, mode2 ###

m = 2
n = 2

set1 = [1,0]
set2 = [0,1]
physical_interferometer = Fourier(m)
part = Partition([Subset(set1), Subset(set2)])

(physical_indexes,  pdf) = compute_probabilities_partition(physical_interferometer, part, n)

print_pdfs(physical_indexes, pdf,n; partition_spans_all_modes = true, physical_events_only = true)

check_photon_conservation(physical_indexes, pdf, n; partition_spans_all_modes = true)
