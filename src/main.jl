using Revise
using BosonSampling
using Permanents
using PrettyTables

# how to store probabilities ? In EventProbability
# array of [counts, proba]
# what is its type?

s1 = Subset([1,1,0,0,0])
s2 = Subset([0,0,1,1,0])
n = 2

part = Partition([s1,s2])
part_occ = PartitionOccupancy(ModeOccupation([2,1]),n,part)

data = MultipleCounts()
