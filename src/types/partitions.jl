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


function partition_from_subset_lengths(subset_lengths)

        """returns a partition from a vector of subset lengths, such as [2,1] gives Partition([[1,1],[1]])"""

    m = sum(subset_lengths)
    subsets = []

    mode = 1

    for subset_length in subset_lengths
        subset_vector = zeros(Int, m)

        for j in 0:subset_length-1
            if mode+j <= m
                subset_vector[mode+j ] = 1
            end
        end

        Subset(subset_vector)
        mode += subset_length
        push!(subsets, Subset(subset_vector))
    end

    Partition(convert(Vector{Subset}, subsets))

end


function equilibrated_partition_vector(m,n_subsets)

    """returns a (most) equilibrated partition possible by euclidian division

    (a problem is that euclidian distribution may give n_subsets or n_subsets+1 if not done like below - here it is the most obvious thing I could think of to get a somewhat equilibrated partition with a constant number of subsets)"""
    q = div(m,n_subsets)
    y = n_subsets
    r = rem(m,n_subsets)

    first_part = [q for i in 1:y]
    first_part[1] += r

    first_part

end

equilibrated_mode_occupation(m,n_subsets) = ModeOccupation(equilibrated_partition_vector(m,n_subsets))

function equilibrated_partition(m,n_subsets)

        """returns a most equilibrated_partition according to the principles of equilibrated_partition_vector"""

        partition_from_subset_lengths(equilibrated_partition_vector(m,n_subsets))

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

                @argcheck counts.m == partition.n_subset "counts do not have as many modes as parition has subsets"
                @argcheck sum(counts.state) <= n "more photons at the output than input"
                if occupies_all_modes(partition)
                        @argcheck sum(counts.state) == n "photons lost in a partition that occupies all modes"
                end

                new(counts, partition, n, partition.subsets[1].m)

        end
end

Base.show(io::IO, part_occ::PartitionOccupancy) = begin

        for (i, count) in enumerate(part_occ.counts.state)
                println(io, count," in ", part_occ.partition.subsets[i])
        end

end

function partition_occupancy_to_partition_size_vector_and_counts(part_occ::PartitionOccupancy)

    part = part_occ.partition

    @argcheck occupies_all_modes(part) "need to have a partition that occupies all modes if using Shsnovitch's formulas in partition_expectation_values.jl"

    partition_size_vector = [part.subsets[i].n for i in 1:length(part.subsets)]
    partition_counts = part_occ.counts.state

    partition_size_vector, partition_counts

end
