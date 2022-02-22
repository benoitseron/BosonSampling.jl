# following Asymptotic Gaussian law for
# noninteracting indistinguishable
# particles in random networks
# Valery S. Shchesnovich

function partition_expectation_values(partition_size_vector, partition_counts)

    """returns the haar averaged probability of photon number count in binned
    outputs for distinguishable and bosonic particles as described in
    Asymptotic Gaussian law for noninteracting indistinguishable particles
    in random networks by Valery S. Shchesnovich (eqs 1,4)"""

    m = sum(partition_size_vector)
    number_partitions = length(partition_size_vector)
    partition_size_ratio = partition_size_vector ./ m # q
    n = sum(partition_counts)
    @test all(partition_counts .>= 0)

    # factorials can quickly get pretty big so this is required
    if any(partition_counts .> 20) || n > 20
        n = big(n)
        partition_counts = big.(partition_counts)
    end

    proba_dist = factorial(n) / vector_factorial(partition_counts) * prod(partition_size_ratio .^ partition_counts)

    proba_bosonic = proba_dist * prod([partition_counts[i] > 1 ? prod([1+l/partition_size_vector[i] for l in 1:partition_counts[i]-1]) : 1 for i in 1:length(partition_size_vector)]) /prod([1+l/m for l in 1:n-1])

    proba_dist, proba_bosonic

end

function subset_expectation_value(subset_size, k,n,m)

    """returns the haar average probability of finding k photons
    inside a subset of binned output modes of size subset_size for the
    distinguishable and bosonic case with n photons and m modes"""

    partition_size_vector = [subset_size, m-subset_size]
    partition_counts = [k,n-k]

    partition_expectation_values(partition_size_vector, partition_counts)

end

function subset_relative_distance_of_averages(subset_size,n,m)

    """returns the distance between 0.5 ∑ₖ|<P_B(k) - P_D(k)>|"""

    @warn "check tvd conventions"
    0.5*abs(sum([abs(subset_expectation_value(subset_size, k,n,m)[1] - subset_expectation_value(subset_size, k,n,m)[2]) for k in 0:n]))


end

function choose_best_average_subset(;m,n, distance = tvd)

    """returns the ideal subset size on average and its TVD,
    to compare with the full bunching test of shsnovitch"""

    function distance_this_subset_size(subset_size)

        proba_dist = [subset_expectation_value(subset_size,k,n,m)[1] for k in 0:n]
        proba_bos = [subset_expectation_value(subset_size,k,n,m)[2] for k in 0:n]

        distance(proba_dist, proba_bos)

    end

    max_distance = 0
    subset_size_max_ratio = nothing

    for subset_size in 1:m-1

        if distance_this_subset_size(subset_size) > max_distance
            max_distance = distance_this_subset_size(subset_size)
            subset_size_max_ratio = subset_size
        end
    end

    subset_size_max_ratio, distance_this_subset_size(subset_size_max_ratio)

end


# 
# ###### averages #######
#
# ### now various plots to see how it evolves ###
#
# n_array = collect(2:20)
#
# m_array = 5 .* n_array
#
# const_density = [choose_best_average_subset(n = n_array[i],m = m_array[i])[2] for i in 1:length(n_array)]
#
# m_array = n_array .^2
#
# no_collision = [choose_best_average_subset(n = n_array[i],m = m_array[i])[2] for i in 1:length(n_array)]
#
# scatter(n_array, const_density)
# scatter!(no_collision)
#
# ### evolution with n ###
#
# n_array = collect(1:40)
#
# m_array = 5 .* n_array
# subset_size_array = map(x-> Int(ceil(x/2)), m_array)
#
# const_density = @. subset_relative_distance_of_averages(subset_size_array, n_array,m_array)
#
# m_array = n_array .^2
# subset_size_array = map(x-> Int(ceil(x/2)), m_array)
#
# no_collision = @. subset_relative_distance_of_averages(subset_size_array, n_array,m_array)
#
# scatter(n_array, const_density)
# scatter!(no_collision)
#
# # this shows that the previous haar averaging was very effective !
#
# ### evolution with partition size ###
#
# n_array = 20
#
# m_array = 5 .* n_array
# subset_size_array = collect(1:m_array-1)
#
# const_density = @. subset_relative_distance_of_averages(subset_size_array, n_array,m_array)
#
# m_array = n_array .^2
# subset_size_array = collect(1:m_array-1)
#
# no_collision = @. subset_relative_distance_of_averages(subset_size_array, n_array,m_array)
#
# scatter(const_density)
# scatter(no_collision)
#
# # this shows that above 1.5 the number of photons, this is pretty much constant
#
# ### evolution with density ###
#
#
# n_array = 20
#
# m_array = collect(n_array:n_array^3)
# subset_size_array = map(x-> Int(ceil(x/2)), m_array)
#
# density_evolution = @. subset_relative_distance_of_averages(subset_size_array, n_array,m_array)
#
# scatter(log10.(density_evolution))
