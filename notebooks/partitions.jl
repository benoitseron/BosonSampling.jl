# this notebook investigates the behaviour of partitions
# and how they can be used for certification

using BosonSampling
using Plots
###### subsets #######

### now various plots to see how it evolves ###

n_array = collect(2:20)

m_array = 5 .* n_array

const_density = [choose_best_average_subset(n = n_array[i],m = m_array[i])[2] for i in 1:length(n_array)]

m_array = n_array .^2

no_collision = [choose_best_average_subset(n = n_array[i],m = m_array[i])[2] for i in 1:length(n_array)]

scatter(n_array, const_density)
scatter!(no_collision)

### evolution with n ###

n_array = collect(1:40)

m_array = 5 .* n_array
subset_size_array = map(x-> Int(ceil(x/2)), m_array)

const_density = @. subset_relative_distance_of_averages(subset_size_array, n_array,m_array)

m_array = n_array .^2
subset_size_array = map(x-> Int(ceil(x/2)), m_array)

no_collision = @. subset_relative_distance_of_averages(subset_size_array, n_array,m_array)

scatter(n_array, const_density)
scatter!(no_collision)

# this shows that the previous haar averaging was very effective !

### evolution with partition size ###

n_array = 20

m_array = 5 .* n_array
subset_size_array = collect(1:m_array-1)

const_density = @. subset_relative_distance_of_averages(subset_size_array, n_array,m_array)

m_array = n_array .^2
subset_size_array = collect(1:m_array-1)

no_collision = @. subset_relative_distance_of_averages(subset_size_array, n_array,m_array)

scatter(const_density)
scatter(no_collision)

# this shows that above 1.5 the number of photons, this is pretty much constant

### evolution with density ###


n_array = 20

m_array = collect(n_array:n_array^3)
subset_size_array = map(x-> Int(ceil(x/2)), m_array)

density_evolution = @. subset_relative_distance_of_averages(subset_size_array, n_array,m_array)

scatter(log10.(density_evolution))

###### multisets ######

m = 12
n = 8
n_subsets_array = collect(2:4)

function plot_tvd_with_n_subsets(m,n,n_subsets_array)

    tvd_n_subsets = [best_partition_size(m = m, n = n, n_subsets = n_subsets)[2] for n_subsets in n_subsets_array]

    plt = scatter(n_subsets_array, tvd_n_subsets)
    xlabel!("n_subsets")
    ylabel!("haar avg tvd")
    title!("m = $m, n = $n")
    savefig("notebooks/images/tvd_with_number_multisets_m=$(m)_n=_$n.png")
    plt
end

plot_tvd_with_n_subsets(m,n,n_subsets_array)
