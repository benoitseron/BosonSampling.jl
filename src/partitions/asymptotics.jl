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

n_detected = 30

bar(sorted_counts[n_detected].proba)