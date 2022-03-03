using Revise
using BosonSampling
using Permanents
using PrettyTables

# how to store probabilities ? In EventProbability
# array of [counts, proba]
# what is its type?

i = first_modes(3,5)

s = Subset([1,1,0,0,0])

Base.show(io::IO, part::Partition) = begin

    println(io, "partition =")
    for s in part.subsets
        println(io, s)
    end
end

s1 = Subset([1,1,0,0,0])
s2 = Subset([0,0,1,1,0])

part = Partition([s1,s2])
