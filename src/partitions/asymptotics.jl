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

# compute partition_size_vector from part
partition_size_vector(part::Partition) = [sum(subset.subset) for subset in part.subsets]


partition_expectation_values_all(part::Partition, )

count = mc.counts[1]

partition_occupancy_array = [PartitionOccupancy(count, part) for count in mc.counts]

# need to only compute the expectation values for even total photon number because of input gaussian states

is_total_photon_number_even(partition_occupancy::PartitionOccupancy) = sum(partition_occupancy.counts.state) % 2 == 0

is_total_photon_number_even.(partition_occupancy_array)


################## needs to be rescaled with overall probability to observe n photons in the whole of the bins

partition_expectation_values_array_dist = [ is_total_photon_number_even(occ) ? partition_expectation_values(occ)[1] : 0 for occ in partition_occupancy_array]

partition_expectation_values_array_bosonic = [ is_total_photon_number_even(occ) ? partition_expectation_values(occ)[2] : 0 for occ in partition_occupancy_array]

bar(partition_expectation_values_array_bosonic, alpha = 0.5)
bar!(partition_expectation_values_array_dist, alpha = 0.5)

