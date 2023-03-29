include("gaussian_partition.jl")


begin

    m = 20
    input_state = GeneralGaussian(m = m, r = 0.40 * ones(m))
    interferometer = RandHaar(m)
    part = equilibrated_partition(m, 2)

    # part = Partition(Subset(first_modes(1, m)))
    # part.subsets[1].subset

    n_max = 50
    mc = compute_probabilities_partition_gaussian(interferometer, part, input_state, n_max)

    # bar(real(pdf))

    sort_samples_total_photon_number_in_partition!(mc)

    @show mc

    bar(mc.proba)

end


### separating events by photon number ###

sorted_counts = sort_by_detected_photons(mc)

n_detected = 14

bar(sorted_counts[n_detected].proba)


### cutoff ###

m = 20
input_state = GeneralGaussian(m = m, r = 1.5 * ones(m))


########### need to be modified - it's for photon pairs !

find_cutoff(input_state; atol = 1e-10)

### average number of photons distribution ###

# m = 20
# input_state = GeneralGaussian(m = m, r = 1.5 * ones(m))

k_range = 0:200

foo(k) = probability_n_photons(k, input_state)

bar(k_range, foo.(k_range))

r_range = 0:0.1:10

plot(r_range, sech.(r_range))
