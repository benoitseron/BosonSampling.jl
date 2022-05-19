using BosonSampling
using Test
using Plots
using JLD

### scattering ###

m = 5
n = 3

interf = RandHaar(m)
i = Input{RandomGramMatrix}(first_modes(n,m))
o = FockDetection(first_modes(n,m))
ev = Event(i,o,interf)

compute_probability!(ev)

ev.proba_params

### bunching ###

m = 5
n = 3

interf = RandHaar(m)
ib = Input{Bosonic}(first_modes(n,m))
ipd = Input{RandomGramMatrix}(first_modes(n,m))
subset_modes = first_modes(n,m)

typeof(subset_modes)

pb = full_bunching_probability(interf, ib, subset_modes)
ppd = full_bunching_probability(interf, ipd, subset_modes)

@test pb/ppd > 1. # this doesn't HAVE TO pass but will pass in nearly all
# cases

### Classical sampling ###

m = 16
n = 3

input = Input{Distinguishable}(first_modes(n,m))
interf = RandHaar(m)
out = classical_sampler(input=input, interf=interf)

## Cliffords sampler ###

m = 16
n = 3
input = Input{Bosonic}(first_modes(n,m))
interf = RandHaar(m)
out = cliffords_sampler(input=input, interf=interf)

## Noisy sampling ###

m = 16
n = 3

x = 0.8 # distinguishability
η = 0.8 # reflectivity

input = Input{OneParameterInterpolation}(first_modes(n,m), x)
interf = RandHaar(m)

out = noisy_sampler(input=input, reflectivity=η, interf=interf)

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
known_sampler = () -> iterate_until_collisionless(() -> classical_sampler(U, m, n)) # gives a classical sampler


samples = metropolis_sampler(;target_pdf = target_pdf, known_pdf = known_pdf , known_sampler = known_sampler , starting_state = starting_state, n_iter = 100)


### Noisy distribution ###

n = 3
m = 6
x = 0.8
η = 0.8

input = Input{OneParameterInterpolation}(first_modes(n,m), x)
interf = RandHaar(m)

output_statistics = noisy_distribution(input=input, reflectivity=η, interf=interf)
p_exact = output_statistics[1]
p_approx = output_statistics[2]
p_sampled = output_statistics[3]

fig_approx = plot(title="approximative computation", xlabel="modes occupation", ylabel="probability");
plot!(fig_approx, p_exact, label="p_exact");
plot!(fig_approx, p_approx, label="p_approx");
fig_samp = plot(title="sampling computation", xlabel="modes occupation", ylabel="probability");
plot!(fig_samp, p_exact, label="p_exact");
plot!(fig_samp, p_sampled, label="p_sampled");

plot(fig_approx, fig_samp, layout=(2,1))

### Theoretical distribution ###

n = 3
m = 6

input = Input{Bosonic}(first_modes(n,m))
interf = RandHaar(m)

output_distribution = theoretical_distribution(input=input, interf=interf)

### Usage Interferometer ###

B = BeamSplitter(1/sqrt(2))
n = 2 # photon number
m = 2 # mode number
proba_bunching = Vector{Float64}(undef, 0)

for x = 0.0:0.01:1.0
    input = Input{OneParameterInterpolation}(first_modes(n,m), x)
    p_theo = theoretical_distribution(input=input, interf=B)

    push!(proba_bunching, p_theo[2]) # store the probabilty to observe one photon in each mode
end
plot(0:0.01:1, proba_bunching, label=nothing, xlabel="distinguishability", ylabel="event probabilty")

### subsets ###

s1 = Subset([1,1,0,0,0])
s2 = Subset([0,0,1,1,0])
s3 = Subset([1,0,1,0,0])

"subsets are not allowed to overlap"

# check_subset_overlap([s1,s2,s3]) will fail

### HOM tests: one mode ###

input_state = Input{Bosonic}(first_modes(n,m))

set1 = [1,0]
physical_interferometer = Fourier(m)
part = Partition([Subset(set1)])

(physical_indexes, pdf) = compute_probabilities_partition(physical_interferometer, part, input_state)


### HOM tests: mode1, mode2 ###

n = 2
m = 2

input_state = Input{Bosonic}(first_modes(n,m))

set1 = [1,0]
set2 = [0,1]
physical_interferometer = Fourier(m)
part = Partition([Subset(set1), Subset(set2)])

(physical_indexes, pdf) = compute_probabilities_partition(physical_interferometer, part, input_state)

print_pdfs(physical_indexes, pdf,n; partition_spans_all_modes = true, physical_events_only = true)

# for a single count

part_occ = PartitionOccupancy(ModeOccupation([1,1]),2,part)

compute_probability_partition_occupancy(physical_interferometer, part_occ, input_state)

# the same using an event

PartitionCount(part_occ)

o = PartitionCount(part_occ)
ev = Event(input_state, o, physical_interferometer)
############ need to change the constructor of Event

get_parametric_type(input_state)
get_parametric_type(o)

length(collect(typeof(o).parameters))
typeof(o)
collect(typeof(input_state).parameters)


### multiset for a random interferometer ###

m = 4
n = 3

inp = Input{Bosonic}(first_modes(n,m))
set1 = zeros(Int,m)
set2 = zeros(Int,m)
set1[1:2] .= 1
set2[3:4] .= 1

physical_interferometer = RandHaar(m)
part = Partition([Subset(set1), Subset(set2)])


(physical_indexes,  pdf) = compute_probabilities_partition(physical_interferometer, part, inp)
fourier_indexes = copy(physical_indexes)


print_pdfs(physical_indexes, pdf, n; physical_events_only = true, partition_spans_all_modes = true)
#print_pdfs(physical_indexes,  probas_fourier, n)

### partitions, subsets ###

n = 2
m = 5

s1 = Subset([1,1,0,0,0])
s2 = Subset([0,0,1,1,0])

part = Partition([s1,s2])
part_occ = PartitionOccupancy(ModeOccupation([2,0]),n,part)

i = Input{Bosonic}(first_modes(n,m))
o = PartitionCount(part_occ)
interf = RandHaar(m)
ev = Event(i,o,interf)

compute_probability!(ev)

### multiple counts probabilities ###


m = 10
n = 3
set1 = zeros(Int,m)
set2 = zeros(Int,m)
set1[1:2] .= 1
set2[3:4] .= 1

interf = RandHaar(m)
part = Partition([Subset(set1), Subset(set2)])

i = Input{Bosonic}(first_modes(n,m))
o = PartitionCountsAll(part)
ev = Event(i,o,interf)

compute_probability!(ev)

### Visualization ###
n = 4
m = 8
input = Input{Bosonic}(ModeOccupation(random_occupancy(n,m)))
interf = RandHaar(m)
output = cliffords_sampler(input=input, interf=interf)
visualize_sampling(input, output)

distinguishability = 0.7
reflectivity = 0.7
input = Input{OneParameterInterpolation}(ModeOccupation(random_occupancy(n,m)), distinguishability)
interf = RandHaar(m)
data_exact = noisy_distribution(input=input, reflectivity=reflectivity, interf=interf, approx=false, samp=false)
data_exact = data_exact[1]

nlist = output_mode_occupation(n,m) # get all the possible output
visualize_proba(input, nlist[200], data_exact)

### Access to datasets ###
data = load("datasets/cliffords_sampler.jld")
output = data["first_modes"]["n=10;m=20;interf=Fourier"]
output = data["random_occupancy"]["n=6;m=36;interf=RandHaar"]
data = load("datasets/noisy_distribution.jld")
p = data["OneParameterInterpolation"]["n=4;m=16;interf=Fourier;reflectivity=0.25;distinguishability_param=0.75"]
p_theo = p[1]
p_approx = p[2]
p_samp = p[3]
