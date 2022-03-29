using Revise
using BosonSampling
using Permanents
using PrettyTables
using ArgCheck

# how to store probabilities ? In EventProbability
# array of [counts, proba]
# what is its type?

s1 = Subset([1,1,0,0,0])
s2 = Subset([0,0,1,1,0])
s3 = Subset([0,0,0,0,1])
n = 3

part = Partition([s1,s2])
part_occ = PartitionOccupancy(ModeOccupation([2,1]),n,part)

occupies_all_modes(part)

part = Partition([s1,s2,s3])
part_occ = PartitionOccupancy(ModeOccupation([2,0,1]),n,part)

partition_expectation_values(part_occ)

m = 7
n = 5
n_subsets = 3
distance = tvd

function best_partition_size(;m,n, n_subsets, distance = tvd)

    """return the ideal partition_size_vector for a given number
    of subsets n_subsets

    for a single subset, n_subsets = 2 as we need a complete partition, occupying all modes"""

    @argcheck n_subsets >= 2 "we need a complete partition, occupying all modes"
    @argcheck n_subsets <= m "more partition bins than output modes"

    part_list = all_mode_configurations(m,n_subsets, only_photon_number_conserving = true)

    remove_trivial_partitions!(part_list)

    max_distance = 0
    part_max_ratio = nothing

    for part in ranked_partition_list(part_list)

        events = all_mode_configurations(m,n_subsets, only_photon_number_conserving = true)

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

best_partition_size(m = m, n = n, n_subsets = n_subsets)


##### eventually make a table of best partition sizes saved in a file ?





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
