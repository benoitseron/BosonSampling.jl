### in this file we compute the probabilities to find
# photons in partitions of the output modes
# a partition is a set of subsets of the output modes

struct Subset
        n::Int
        m::Int
        state::Vector{Int}
        Subset(state) = isa_subset(state) ? new(sum(state), length(state), state) : error("invalid subset")
end

struct Partition
        subsets::Vector{Subset}
        n_subset::Int
        Partition(subsets) = new(subsets, n_subset)
end

struct PartitionOccupancy
        counts::ModeOccupation
        partition::Partition
        function PartitionOccupancy(counts, partition)
                if counts.m == partition.n_subset
                        new(counts, partition)
                else
                        error("counts do not have as many modes as parition has subsets")
                end
        end
end
