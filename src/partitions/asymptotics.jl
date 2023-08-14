include("gaussian_partition.jl")

### binning demo ###

# overall bins:

begin

    m = 20
    input_state = GeneralGaussian(m = m, r = 0.4 * ones(m))
    interferometer = RandHaar(m)
    part = equilibrated_partition(m, 2)

    # part = Partition(Subset(first_modes(1, m)))
    # part.subsets[1].subset

    n_max = 10
    mc = compute_probabilities_partition_gaussian(interferometer, part, input_state, n_max)

    # bar(real(pdf))

    sort_samples_total_photon_number_in_partition!(mc)

    @show mc

    bar(mc.proba)

end

bar(mc.proba[1:100])

mc.counts

# zooming

begin

    m = 20
    input_state = GeneralGaussian(m = m, r = 0.4 * ones(m))
    interferometer = RandHaar(m)
    part = equilibrated_partition(m, 2)

    # part = Partition(Subset(first_modes(1, m)))
    # part.subsets[1].subset

    n_max = 40
    mc = compute_probabilities_partition_gaussian(interferometer, part, input_state, n_max)

    # bar(real(pdf))

    sort_samples_total_photon_number_in_partition!(mc)

    @show mc

    bar(mc.proba)

end


### separating events by photon number ###

sorted_counts = sort_by_detected_photons(mc)

n_detected = 24

bar(sorted_counts[n_detected].proba)


### plotting the average number of photons distribution ###


n_max = 10

if n_max % 2 == 0
    @warn "n_max must be odd for FFT purposes, converting"
    n_max = n_max +1
end

counts = all_mode_configurations(n_max,part.n_subset, only_photon_number_conserving = false)



foo()

"""
    function partition_expectation_values_gaussian(input_state::GeneralGaussian, n_max::Int, part::Partition)

Computes the asymptotic law of the expectation values of a Haar averaged partition for Gaussian input states .
"""
function partition_expectation_values_gaussian(input_state::GeneralGaussian, n_max::Int, part::Partition)


    #generating the list of indexes for counts 
    if n_max % 2 == 0
        @warn "n_max must be odd for FFT purposes, converting"
        n_max = n_max +1
    end

    counts = all_mode_configurations(n_max,part.n_subset, only_photon_number_conserving = false)

    partition_occupancy_array = [PartitionOccupancy(count, part) for count in counts]

    # need to only compute the expectation values for even total photon number because of input gaussian states

    is_total_photon_number_even(partition_occupancy::PartitionOccupancy) = sum(partition_occupancy.counts.state) % 2 == 0

    partition_expectation_values_array_dist = [is_total_photon_number_even(occ) ? probability_n_photons(occ.n, input_state) * partition_expectation_values(occ)[1] : 0 for occ in partition_occupancy_array]

    partition_expectation_values_array_bosonic = [is_total_photon_number_even(occ) ? probability_n_photons(occ.n, input_state) * partition_expectation_values(occ)[2] : 0 for occ in partition_occupancy_array]

    partition_expectation_values_array_bosonic, partition_expectation_values_array_dist

end


partition_expectation_values_array_bosonic, partition_expectation_values_array_dist = partition_expectation_values_gaussian(input_state, n_max, part)

bar(partition_expectation_values_array_bosonic, alpha = 0.3, label = L"Indistinguishable")
bar!(partition_expectation_values_array_dist, alpha = 0.3, label = L"Distinguishable")

##### THIS MIGHT EXPLAIN THE MISTAKE WITH THE CHICAGO PEOPLE #####

# what changes compared to the previous commit is that I use 

#  #generating the list of indexes for counts 
#  if n_max % 2 == 0
#     @warn "n_max must be odd for FFT purposes, converting"
#     n_max = n_max +1
# end

# counts = all_mode_configurations(n_max,part.n_subset, only_photon_number_conserving = false)

# and these look more like the chicago people's code, and I think it's simply because the list was no re-ordered by total photon number