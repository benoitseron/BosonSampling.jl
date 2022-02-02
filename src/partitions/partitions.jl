### in this file we compute the probabilities to find
# photons in partitions of the output modes
# a partition is a set of subsets of the output modes

struct Subset
        # basically a mode occupation list with at most one count per mode
        n::Int
        m::Int
        subset::Vector{Int}
        Subset(state) = isa_subset(state) ? new(sum(state), length(state), state) : error("invalid subset")
end

struct Partition
        subsets::Vector{Subset}
        n_subset::Int
        Partition(subsets) = new(subsets, length(subsets))
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

struct PartitionCountingInterferometer <: Interferometer
        partition_occupancy::PartitionOccupancy
        physical_interferometer::Interferometer
        fourier_indexes::Vector{Complex}
        virtual_interferometer::Interferometer
        # the phases are 2pi/(n+1) * fourier_indexes

        function PartitionCountingInterferometer(partition_occupancy, physical_interferometer, fourier_indexes)

                """outputs the virtual interferometer that gives the fourier phase
                x(fourier_indexes), of which the multidimensional inverse
                fourier transform outputs the the probability to get partition_occupancy.counts
                photons in partition_occupancy.partition"""

                virtual_interferometer_matrix = copy(physical_interferometer.U)
                for i in length(partition_occupancy.partition.subsets)
                        virtual_interferometer_matrix *= Diagonal(@. exp(2*pi*1im/(partition_occupancy.n+1) * fourier_indexes))
                end
                virtual_interferometer_matrix *= U'

                new(partition_occupancy, physical_interferometer, virtual_interferometer_matrix, fourier_indexes, virtual_interferometer_matrix)
        end

end

function all_mode_configurations(n,n_subset; only_photon_number_conserving = false)

        """generates all possible output modes configurations for n photons
        in n_subset outputs

        does not take into account photon number conservation by default"""

        array = []
        for i in 1:(n+1)^(n_subset)

                this_vector = digits(i-1, base = n+1, pad = n_subset)

                if only_photon_number_conserving
                        if sum(this_vector) == n
                                push!(array,this_vector)
                        end
                else
                        push!(array,this_vector)
                end

        end
        array

end

m = 6
n = 4
set1 = zeros(Int,m)
set2 = zeros(Int,m)
set1[1:2] .= 1
set2[3:4] .= 1

physical_interferometer = RandHaar(m)
part = Partition([Subset(set1), Subset(set2)])
partition_occupancy = PartitionOccupancy(ModeOccupation([1,2]), n, part)


fourier_indexes = all_mode_configurations(n,part.n_subset)
probas_fourier = Array{ComplexF64}(undef, length(fourier_indexes))
virtual_interferometer_matrix = similar(physical_interferometer.U)
#
# for (i, fourier_index) in enumerate(fourier_indexes)
#         probas_fourier[i] = permanent(PartitionCountingInterferometer(partition_occupancy, physical_interferometer, fourier_index).virtual_interferometer.U)
# end
#
# virtual_interferometer_matrix = copy(physical_interferometer.U)
# v = 1
#
# this_diag = [one(eltype(virtual_interferometer_matrix)) for i in 1:m]
# this_diag *= subset .* this_phase
#

partition_occupancy.partition.subsets

for (index_fourier_array, fourier_index) in enumerate(fourier_indexes)

        # for each fourier index, we recompute the virtual interferometer
        virtual_interferometer_matrix  = physical_interferometer.U'
        diag = [one(eltype(virtual_interferometer_matrix)) for i in 1:m]

        @show fourier_index
        for (i,fourier_element) in enumerate(fourier_index)

                this_phase = exp(2*pi*1im/(partition_occupancy.n+1) * fourier_element)
                @show this_phase
                @show partition_occupancy.partition.subsets[i].subset
                for j in 1:length(diag)
                        @show j
                        if partition_occupancy.partition.subsets[i].subset[j] == 1
                                diag[j] *= this_phase
                                @show diag[j]
                        end
                end

        end
        @show diag
        virtual_interferometer_matrix *= Diagonal(diag)
        virtual_interferometer_matrix *= physical_interferometer.U
        #@show virtual_interferometer_matrix

        probas_fourier[index_fourier_array] = permanent(virtual_interferometer_matrix)
end

probas_fourier
new(partition_occupancy, physical_interferometer, virtual_interferometer_matrix, fourier_indexes, virtual_interferometer_matrix)

exp.([1,2])

v
