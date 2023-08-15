include("gaussian_partition.jl")

using BosonSampling
using Plots
using ProgressMeter
using Distributions
using Random
using Test
using ArgCheck
using StatsBase
using ColorSchemes
using Interpolations
using Dierckx
using LinearAlgebra
using PrettyTables
using LaTeXStrings
using JLD
using AutoHashEquals
using LinearRegression

using DataStructures

color_map = ColorSchemes.rainbow

### interpolation functions for plots ###


"""

    interpolate_data(x_data, y_data, n_points = 1000)

Provides cubic spline interpolation of data for plotting.
"""
function interpolate_data(x_data, y_data, n_points = 1000)
    
    # converting x_data to numerical values
    if !isa(x_data[1], Number)

        x_data = 1:length(x_data)

    end

    x_spl = range(1, stop=length(x_data), length=n_points)
    spl = Spline1D(x_data,y_data)
    y_spl = spl(x_spl)

    return x_spl, y_spl

end


### Haar averaged binning demo ###

n_iter = 10

m = 20
input_state = GeneralGaussian(m = m, r = 0.4 * ones(m))
part = equilibrated_partition(m, 2)

n_max = 20

begin
    
    mc_array = []

    @showprogress for i in 1:n_iter

        interferometer = RandHaar(m)
        mc_iter = compute_probabilities_partition_gaussian(interferometer, part, input_state, n_max)

        # bar(real(pdf))

        sort_samples_total_photon_number_in_partition!(mc_iter)
        push!(mc_array, mc_iter)

    end

end

averaged_proba = vec(mean(hcat([mc_array[i].proba for i in 1:n_iter]), dims = 1))[1]
deviation = vec(sqrt.(var(hcat([mc_array[i].proba for i in 1:n_iter]...), dims = 2)))

mc = MultipleCounts(mc_array[1].counts, averaged_proba)
mc_deviation = MultipleCounts(mc_array[1].counts, deviation)

mc_partition_expectation_values_array_bosonic, mc_partition_expectation_values_array_dist = partition_expectation_values_gaussian(input_state, n_max, part)

# bar(mc_partition_expectation_values_array_bosonic.proba, alpha = 0.5, label = L"Indistinguishable")
# bar!(mc_partition_expectation_values_array_dist.proba, alpha = 0.5, label = L"Distinguishable")

# ### comparison with numerics ###

# mc, probas_fourier = compute_probabilities_partition_gaussian_chicago(interferometer, part, input_state, n_max, give_debug_info = true)

# sort_samples_total_photon_number_in_partition!(mc)

# range_values = 1:50

# bar(real.(mc.proba[range_values]), alpha = 0.5, label = L"Numerics")
# bar!(mc_partition_expectation_values_array_bosonic.proba[range_values], alpha = 0.5, label = L"Indistinguishable")
# bar!(mc_partition_expectation_values_array_dist.proba[range_values], alpha = 0.5, label = L"Distinguishable")

# ### separating events by photon number ###

# sorted_counts_numerics = sort_by_detected_photons(mc)
# sorted_counts_bosonic = sort_by_detected_photons(mc_partition_expectation_values_array_bosonic)
# sorted_counts_dist = sort_by_detected_photons(mc_partition_expectation_values_array_dist)

# n_detected = 10

# scatter(sorted_counts_numerics[n_detected].proba, alpha = 0.5, label = L"Numerics")
# scatter!(sorted_counts_bosonic[n_detected].proba, alpha = 0.5, label = L"Indistinguishable")
# scatter!(sorted_counts_dist[n_detected].proba, alpha = 0.5, label = L"Distinguishable")

# ###### to plot a mc ######


# dpi_value = 600
# fig_width = 1200  # in pixels
# fig_height = 500  # in pixels


# range_values = 1:80
# mc_reduced = MultipleCounts(mc.counts[range_values], mc.proba[range_values])

# x_data = mc_reduced.counts
# y_data = mc_reduced.proba

# function plot_string_repr(i::ModeOccupation)
#     return string(i.state)
# end

# x_labels = map(plot_string_repr, x_data)

# plot(y_data, xticks=(1:length(x_labels), x_labels), xrotation=60, legend=false, xtickfontsize=8, dpi=dpi_value, size=(fig_width, fig_height))



### large plot ###

c1 = :cyan
c2 = :orange
c3 = :purple


marker = :circle

dpi_value = 600
fig_width = 1500  # in pixels
fig_height = 900  # in pixels

range_values = 1:45

mc_reduced = MultipleCounts(mc.counts[range_values], mc.proba[range_values])

x_data = mc_reduced.counts
y_data = mc_reduced.proba
y_err_data = mc_deviation.proba[range_values]

# interpolations

function plot_string_repr(i::ModeOccupation)
    return string(i.state)
end

x_labels = map(plot_string_repr, x_data)

bar(y_data, yerr = y_err_data, xticks=(1:length(x_labels), x_labels), xrotation=60, xtickfontsize=8, ytickfontsize=12, dpi=dpi_value, size=(fig_width, fig_height), alpha = 0.5, label = L"Numerics", c = c1)

scatter!(mc_partition_expectation_values_array_bosonic.proba[range_values], alpha = 1, label = L"Indistinguishable", c = c2, m = marker)

# scatter!(mc_partition_expectation_values_array_dist.proba[range_values], alpha = 1, label = L"Distinguishable", c = c3, m = marker)

### zoomed plot ###

sorted_counts_numerics = sort_by_detected_photons(mc)
sorted_counts_numerics_deviation = sort_by_detected_photons(mc_deviation)
sorted_counts_bosonic = sort_by_detected_photons(mc_partition_expectation_values_array_bosonic)
sorted_counts_dist = sort_by_detected_photons(mc_partition_expectation_values_array_dist)

n_detected = 18

scatter(sorted_counts_numerics[n_detected].proba, alpha = 0.5, label = L"Numerics")
scatter!(sorted_counts_bosonic[n_detected].proba, alpha = 0.5, label = L"Indistinguishable")
scatter!(sorted_counts_dist[n_detected].proba, alpha = 0.5, label = L"Distinguishable")

mc_reduced = sorted_counts_numerics[n_detected]
mc_reduced_deviation = sorted_counts_numerics_deviation[n_detected]
mc_reduced_bosonic = sorted_counts_bosonic[n_detected]
mc_reduced_dist = sorted_counts_dist[n_detected]

x_data = mc_reduced.counts
y_data = mc_reduced.proba
y_err_data = mc_reduced_deviation.proba

function plot_string_repr(i::ModeOccupation)
    return string(i.state)
end

x_labels = map(plot_string_repr, x_data)

# Define margins: [left, bottom, right, top]

bar(y_data, yerr = y_err_data, xticks=(1:length(x_labels), x_labels), xrotation=60, xtickfontsize=8, ytickfontsize=12, dpi=dpi_value, size=(fig_width, fig_height), alpha = 0.5, label = L"Numerics", c = c1)

scatter!(mc_reduced_bosonic.proba, alpha = 1, label = L"Indistinguishable" , c = c2, m = marker)
interpolated_bosonic = interpolate_data(1:length(mc_reduced_bosonic.proba), mc_reduced_bosonic.proba)
plot!(interpolated_bosonic[1], interpolated_bosonic[2], alpha = 1, label = L"Indistinguishable" , c = c2)

# the same for distinguishable

scatter!(mc_reduced_dist.proba, alpha = 1, label = L"Distinguishable", c = c3, m = marker)
interpolated_dist = interpolate_data(1:length(mc_reduced_dist.proba), mc_reduced_dist.proba)
plot!(interpolated_dist[1], interpolated_dist[2], alpha = 1, label = L"Distinguishable", c = c3)


xlabel!(L"configuration")
ylabel!(L"p")







