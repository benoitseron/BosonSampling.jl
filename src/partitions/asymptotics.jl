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


dpi_value = 600
fig_width = 1200  # in pixels
fig_height = 500  # in pixels


range = 1:80
mc_reduced = MultipleCounts(mc.counts[range], mc.proba[range])

x_data = mc_reduced.counts
y_data = mc_reduced.proba

function plot_string_repr(i::ModeOccupation)
    return string(i.state)
end

x_labels = map(plot_string_repr, x_data)

plot(y_data, xticks=(1:length(x_labels), x_labels), xrotation=60, legend=false, xtickfontsize=8, dpi=dpi_value, size=(fig_width, fig_height))

### large plot ###

dpi_value = 600
fig_width = 1500  # in pixels
fig_height = 800  # in pixels

range = 1:45

mc_reduced = MultipleCounts(mc.counts[range], mc.proba[range])

x_data = mc_reduced.counts
y_data = mc_reduced.proba

function plot_string_repr(i::ModeOccupation)
    return string(i.state)
end

x_labels = map(plot_string_repr, x_data)

# Define margins: [left, bottom, right, top]

scatter(y_data, xticks=(1:length(x_labels), x_labels), xrotation=60, legend=false, xtickfontsize=8, ytickfontsize=12, dpi=dpi_value, size=(fig_width, fig_height), alpha = 0.5, label = L"Numerics")
bar!(mc_partition_expectation_values_array_bosonic.proba[range], alpha = 0.5, label = L"Indistinguishable")
bar!(mc_partition_expectation_values_array_dist.proba[range], alpha = 0.5, label = L"Distinguishable")

