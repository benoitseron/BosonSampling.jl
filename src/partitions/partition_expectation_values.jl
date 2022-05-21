# following Asymptotic Gaussian law for
# noninteracting indistinguishable
# particles in random networks
# Valery S. Shchesnovich


"""
    partition_expectation_values(partition_size_vector, partition_counts)
    partition_expectation_values(part_occ::PartitionOccupancy)

Return the Haar averaged probability of photon number count in binned outputs
for [`Distinguishable`](@ref) and [`Bosonic`](@ref) particles.

!!! note "Reference"
        [https://www.nature.com/articles/s41598-017-00044-8.pdf](https://www.nature.com/articles/s41598-017-00044-8.pdf)
"""
function partition_expectation_values(partition_size_vector, partition_counts)

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

partition_expectation_values(part_occ::PartitionOccupancy) = partition_expectation_values(partition_occupancy_to_partition_size_vector_and_counts(part_occ)...)

"""
    subset_expectation_value(subset_size, k, n, m)

Return the Haar averaged probability to find `k` from `n` photons inside a subset of
binned output modes of length `size_subset` among `m` modes for [`Distinguishable`](@ref)
and [`Bosonic`](@ref) cases.
"""
function subset_expectation_value(subset_size, k,n,m)

    partition_size_vector = [subset_size, m-subset_size]
    partition_counts = [k,n-k]

    partition_expectation_values(partition_size_vector, partition_counts)

end

"""
    subset_relative_distance_of_averages(subset_size, n, m)
"""
function subset_relative_distance_of_averages(subset_size,n,m)

    @warn "check tvd conventions"
    0.5*abs(sum([abs(subset_expectation_value(subset_size, k,n,m)[1] - subset_expectation_value(subset_size, k,n,m)[2]) for k in 0:n]))

end

"""
    choose_best_average_subset(;m, n, distance=tvd)

Return the ideal subset size on average and its total variance distance.
"""
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


"""
    best_partition_size(;m, n, n_subsets, distance=tvd)

Return the ideal `partition_size_vector` for a given number of subsets `n_subsets`
and the Haar averaged TVD in second parameter.

!!! note
    For a single subset, `n_subsets`=2 as we need a complete partition, occupying all modes.
"""
function best_partition_size(;m,n, n_subsets, distance = tvd)

    """return the ideal partition_size_vector for a given number
    of subsets n_subsets as well as the haar averaged tvd in second parameter

    for a single subset, n_subsets = 2 as we need a complete partition, occupying all modes"""

    @argcheck n_subsets >= 2 "we need a complete partition, occupying all modes"
    @argcheck n_subsets <= m "more partition bins than output modes"

    part_list = all_mode_configurations(m,n_subsets, only_photon_number_conserving = true)

    remove_trivial_partitions!(part_list)

    max_distance = 0
    part_max_ratio = nothing

    for part in ranked_partition_list(part_list)

        events = all_mode_configurations(n,n_subsets, only_photon_number_conserving = true)

        pdf = [partition_expectation_values(part, event) for event in events]

        pdf_dist = hcat(collect.(pdf)...)[1,:]
        pdf_bos = hcat(collect.(pdf)...)[2,:]

        this_distance = distance(pdf_bos,pdf_dist)

        println(part, " : ", this_distance)

        if this_distance > max_distance
            max_distance = this_distance
            part_max_ratio = part
        end
    end

    part_max_ratio, max_distance

end
