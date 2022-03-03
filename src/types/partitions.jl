### Partitions ###

struct ModeOccupation
    n::Int
    m::Int
    state::Vector{Int}
    ModeOccupation(state) = all(state[:] .>= 0) ? new(sum(state), length(state), state) : error("negative photon counts")
end

Base.show(io::IO, i::ModeOccupation) = print(io, "state = ", i.state)

at_most_one_photon_per_bin(state) = all(state[:] .<= 1)
at_most_one_photon_per_bin(r::ModeOccupation) = at_most_one_photon_per_bin(r.state)

isa_subset(subset_modes::Vector{Int}) = (at_most_one_photon_per_bin(subset_modes) && sum(subset_modes) != 0)
isa_subset(subset_modes::ModeOccupation) = isa_subset(subset_modes.state)

first_modes(n::Int,m::Int) = n<=m ? ModeOccupation([i <= n ? 1 : 0 for i in 1:m]) : error("n>m")


struct Subset
        # basically a mode occupation list with at most one count per mode
        n::Int
        m::Int
        subset::Vector{Int}
        function Subset(state)

                isa_subset(state) ? new(sum(state), length(state), state) : error("invalid subset")
        end
end

Base.show(io::IO, s::Subset) = print(io, "subset = ", convert(Vector{Int},occupancy_vector_to_partition(s.subset)))


function check_disjoint_subsets(s1::Subset, s2::Subset)
        @argcheck s1.m == s2.m "subsets do not have the same dimension"
        @argcheck all(s1.subset .* s2.subset .== 0) "subsets overlap"
end

function check_subset_overlap(subsets::Vector{Subset})

        for (i,subset_1) in enumerate(subsets)
                for (j,subset_2) in enumerate(subsets)
                        if i>j
                                check_disjoint_subsets(subset_1, subset_2)
                        end
                end
        end

end

struct Partition
        subsets::Vector{Subset}
        n_subset::Int
        m::Int
        function Partition(subsets)
                check_subset_overlap(subsets)
                new(subsets, length(subsets), subsets[1].m)
        end
end

Base.show(io::IO, part::Partition) = begin

    println(io, "partition =")
    for s in part.subsets
        println(io, s)
    end
end


function occupies_all_modes(part::Partition)

        """checks if a partition occupies all m modes"""

        occupied_modes = zeros(Int, part.m)
        for s in part.subsets
                occupied_modes .+= s.subset
        end

        all(occupied_modes .== 1)

end

struct PartitionOccupancy
        counts::ModeOccupation
        partition::Partition
        n::Int
        m::Int
        function PartitionOccupancy(counts::ModeOccupation, n::Int, partition::Partition)
                if counts.m == partition.n_subset
                        new(counts, partition, n, partition.subsets[1].m)
                else
                        error("counts do not have as many modes as parition has subsets")
                end
        end
end

Base.show(io::IO, part_occ::PartitionOccupancy) = begin

    println(io, "counts = ", part_occ.counts.state, " in")
    println(io,part_occ.partition)

end
