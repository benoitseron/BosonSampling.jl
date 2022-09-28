"""
    ModeOccupation(state)

A list of the size of the number of modes `m`, with entry `j` of `state` being the number of photons in mode `j`. See also [`ModeList`](@ref).

    fields:
         - n::Int
         - m::Int
         - state::Vector{Int}
"""
@auto_hash_equals struct ModeOccupation
    n::Int
    m::Int
    state::Vector{Int}
    ModeOccupation(state) = all(state[:] .>= 0) ? new(sum(state), length(state), state) : error("negative photon counts")
end

Base.show(io::IO, i::ModeOccupation) = print(io, "state = ", i.state)

"""

    :+(s1::ModeOccupation, s2::ModeOccupation)
    :+(s1::ModeOccupation, s2::Vector{Int})
    :+(s2::Vector{Int}, s1::ModeOccupation)

Adds two mode occupations, for instance
s1 = ModeOccupation([0,1])
s2 = ModeOccupation([1,0])

(s1+s2).state == [1,1]

Also works with just a vector and a mode occupation.
"""
Base.:+(s1::ModeOccupation, s2::ModeOccupation) = begin
    return ModeOccupation(s1.state + s2.state)
end

Base.:+(s1::ModeOccupation, s2::Vector{Int}) = begin

        @argcheck size(s1.state) == size(s2) "incompatible sizes"
    return ModeOccupation(s1.state + s2)
end

Base.:+(s2::Vector{Int}, s1::ModeOccupation) = begin
    return s1 + s2
end

"""
        Base.zeros(mo::ModeOccupation)

Returns a `ModeOccupation` similar to the input but with a state made of zeros.
"""
function Base.zeros(mo::ModeOccupation)

    physical_state = mo.state
    state = zeros(eltype(physical_state), size(physical_state))
    ModeOccupation(state)

end

"""
        Base.cat(s1::ModeOccupation, s2::ModeOccupation)

Concatenates two `ModeOccupation`.
"""
function Base.cat(s1::ModeOccupation, s2::ModeOccupation)

    ModeOccupation(vcat(s1.state, s2.state))

end

"""
    ModeList(state)
    ModeList(state,m)

Contrasting to [`ModeOccupation`](@ref) this list is of size `n`, the number of photons. Entry `j` is the index of the mode occupied by photon `j`.

This can also be used just to select modes for instance.

See also [`ModeOccupation`](@ref).

    fields:
        - n::Int
        - m::Union{Int, Nothing}
        - modes::Vector{Int}
"""
@auto_hash_equals struct ModeList
    n::Int
    m::Union{Int, Nothing}
    modes::Vector{Int}

    ModeList(modes::Vector{Int}) = ModeList(modes, nothing)

    # all(modes[:] .>= 1) ? new(length(modes), nothing, modes) : error("modes start at one")

        function ModeList(modes::Vector{Int}, m)

                if all(modes[:] .>= 1) && (m == nothing ? true : all(modes[:] .<= m))
                    new(length(modes), m, modes)
                else
                    error("incoherent or invalid mode inputs")
                end
        end

        ModeList(mode::Int, m = nothing) = ModeList([mode],m)

end

"""
        is_compatible(target_modes_in::ModeList, target_modes_out::ModeList)

Checks compatibility of `ModeList`s.
"""
function is_compatible(target_modes_in::ModeList, target_modes_out::ModeList)

        if target_modes_in == target_modes_out
                return true
        else
                @argcheck target_modes_in.n == target_modes_out.n
                @argcheck target_modes_in.m == target_modes_out.m
                true
        end

end

function Base.convert(::Type{ModeOccupation}, ml::ModeList)

        if ml.m == nothing
                error("need to give m")
        else
                state = zeros(Int, ml.m)

                for mode in ml.modes
                        state[mode] += 1
                end

                return ModeOccupation(state)
        end
end

at_most_one_photon_per_bin(state) = all(state[:] .<= 1)
at_most_one_photon_per_bin(r::ModeOccupation) = at_most_one_photon_per_bin(r.state)

isa_subset(subset_modes::Vector{Int}) = (at_most_one_photon_per_bin(subset_modes) && sum(subset_modes) != 0)
isa_subset(subset_modes::ModeOccupation) = isa_subset(subset_modes.state)

"""
    first_modes(n::Int, m::Int)

Create a [`ModeOccupation`](@ref) with `n` photons in the first sites of `m` modes.
"""
first_modes(n::Int,m::Int) = n<=m ? ModeOccupation([i <= n ? 1 : 0 for i in 1:m]) : error("n>m")

first_modes_array(n::Int,m::Int) = first_modes(n,m).state

"""
    last_modes(n::Int, m::Int)

Create a [`ModeOccupation`](@ref) with `n` photons in the last sites of `m` modes.
"""
last_modes(n::Int,m::Int) = n<=m ? ModeOccupation([i > m-n ? 1 : 0 for i in 1:m]) : error("n>m")

last_modes_array(n::Int,m::Int) = last_modes(n,m).state


"""
    Subset(state::Vector{Int})

Create a mode occupation list with at most one count per mode.

    Fields:
         - n::Int
         - m::Int
         - subset::Vector{Int}
"""
@auto_hash_equals struct Subset
        # basically a mode occupation list with at most one count per mode
        n::Int
        m::Int
        subset::Vector{Int}
        function Subset(state)

                isa_subset(state) ? new(sum(state), length(state), state) : error("invalid subset")
        end
        function Subset(modeocc::ModeOccupation)

                state = modeocc.state
                isa_subset(state) ? new(sum(state), length(state), state) : error("invalid subset")
        end

        function Subset(ml::ModeList)
                Subset(convert(ModeOccupation, ml))
        end
end

Base.show(io::IO, s::Subset) = print(io, "subset = ", convert(Vector{Int},occupancy_vector_to_partition(s.subset)))


function check_disjoint_subsets(s1::Subset, s2::Subset)
        @argcheck s1.m == s2.m "subsets do not have the same dimension"
        @argcheck all(s1.subset .* s2.subset .== 0) "subsets overlap"
end

function check_subset_overlap(subsets::Vector{Subset})

        if length(subsets) == 1
                return false
        end

        for (i,subset_1) in enumerate(subsets)
                for (j,subset_2) in enumerate(subsets)
                        if i>j
                                check_disjoint_subsets(subset_1, subset_2)
                        end
                end
        end

end

function check_subset_overlap(subset::Subset)
        nothing
end

"""
    Partition(subsets::Vector{Subset})

Create a partition from multiple [`Subset`](@ref).
"""
@auto_hash_equals struct Partition
        subsets::Vector{Subset}
        n_subset::Int
        m::Int
        function Partition(subsets)
                check_subset_overlap(subsets)
                new(subsets, length(subsets), subsets[1].m)
        end

        function Partition(subset::Subset)
                Partition([subset])
        end
end

Base.show(io::IO, part::Partition) = begin

    println(io, "partition =")
    for s in part.subsets
        println(io, s)
    end
end

"""
    partition_from_subset_lengths(subset_lengths)

Return a partition from a vector of subset lengths.
"""
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

"""
    equilibrated_partition_vector(m,n_subsets)

Returns a (most) equilibrated partition possible by euclidian division.

(a problem is that euclidian distribution may give n_subsets or n_subsets+1 if not done like below - here it is the most obvious thing I could think of to get a somewhat equilibrated partition with a constant number of subsets)
"""
function equilibrated_partition_vector(m,n_subsets)


    q = div(m,n_subsets)
    y = n_subsets
    r = rem(m,n_subsets)

    first_part = [q for i in 1:y]
    first_part[1] += r

    first_part

end

equilibrated_mode_occupation(m,n_subsets) = ModeOccupation(equilibrated_partition_vector(m,n_subsets))

"""
    equilibrated_partition(m,n_subsets)

Returns a most equilibrated_partition according to the principles of [`equilibrated_partition_vector`](@ref).
"""
function equilibrated_partition(m,n_subsets)

        partition_from_subset_lengths(equilibrated_partition_vector(m,n_subsets))

end

"""
    occupies_all_modes(part::Partition)

Check wether a partition occupies all modes or not.
"""
function occupies_all_modes(part::Partition)

        """checks if a partition occupies all m modes"""

        occupied_modes = zeros(Int, part.m)
        for s in part.subsets
                occupied_modes .+= s.subset
        end

        all(occupied_modes .== 1)

end

"""
    PartitionOccupancy(counts::ModeOccupation, n::Int, partition::Partition)

   Fields:
        - counts::ModeOccupation
        - partition::Partition
        - n::Int
        - m::Int
"""
@auto_hash_equals struct PartitionOccupancy
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

remove_last_subset(part::Partition) = Partition(part.subsets[1:end-1])
