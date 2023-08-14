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

    n_max = 30
    mc = compute_probabilities_partition_gaussian(interferometer, part, input_state, n_max)

    # bar(real(pdf))

    sort_samples_total_photon_number_in_partition!(mc)

    @show mc

    bar(mc.proba)

end


mc_partition_expectation_values_array_bosonic, mc_partition_expectation_values_array_dist = partition_expectation_values_gaussian(input_state, n_max, part)

bar(mc_partition_expectation_values_array_bosonic.proba, alpha = 0.5, label = L"Indistinguishable")
bar!(mc_partition_expectation_values_array_dist.proba, alpha = 0.5, label = L"Distinguishable")

### comparison with numerics ###

mc, probas_fourier = compute_probabilities_partition_gaussian_chicago(interferometer, part, input_state, n_max, give_debug_info = true)

sort_samples_total_photon_number_in_partition!(mc)

range = 1:50

bar(real.(mc.proba[range]), alpha = 0.5, label = L"Numerics")
bar!(mc_partition_expectation_values_array_bosonic.proba[range], alpha = 0.5, label = L"Indistinguishable")
bar!(mc_partition_expectation_values_array_dist.proba[range], alpha = 0.5, label = L"Distinguishable")

### separating events by photon number ###

sorted_counts_numerics = sort_by_detected_photons(mc)
sorted_counts_bosonic = sort_by_detected_photons(mc_partition_expectation_values_array_bosonic)
sorted_counts_dist = sort_by_detected_photons(mc_partition_expectation_values_array_dist)

n_detected = 10

scatter(sorted_counts_numerics[n_detected].proba, alpha = 0.5, label = L"Numerics")
scatter!(sorted_counts_bosonic[n_detected].proba, alpha = 0.5, label = L"Indistinguishable")
scatter!(sorted_counts_dist[n_detected].proba, alpha = 0.5, label = L"Distinguishable")

###### to plot a mc ######

x_data = ["string1", "string2", "string3"]
y_data = [1, 2, 3]

plot(y_data, xticks=(1:length(x_data), x_data), xrotation=45, legend=false, xguidefontsize=8)
