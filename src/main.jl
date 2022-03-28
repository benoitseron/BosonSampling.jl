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

partition_expectation_values(part_data...)
partition_expectation_values(part_occ)
