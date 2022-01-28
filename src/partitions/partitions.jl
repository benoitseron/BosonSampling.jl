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
        function PartitionOccupancy(counts, partition)
                if counts.m == partition.n_subset
                        new(counts, partition, counts.n, subset.m)
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
                for i in length(partition_occupancy.subsets)
                        virtual_interferometer_matrix *= Diagonal(@. exp(2*pi*1im/(partition_occupancy.n+1) * fourier_indexes))
                end
                virtual_interferometer_matrix *= U'

                new(partition_occupancy, physical_interferometer, virtual_interferometer_matrix, fourier_indexes, virtual_interferometer_matrix)
        end

end

m = 6
n = 4
set1 = zeros(Int,m)
set2 = zeros(Int,m)
set1[1:2] .= 1
set2[3:4] .= 1

part = Partition([Subset(set1), Subset(set2)])

p = 5


function all_mode_configuration(n,n_subset; only_photon_number_conserving = false)

        """prints all possible output modes configurations for n photons
        in n_subset outputs

        does not take into account photon number conservation by default"""
        for i in 1:(n+1)^(n_subset)

                this_vector = digits(i-1, base = n+1, pad = n_subset)

                if only_photon_number_conserving
                        if sum(this_vector) == n
                                println(this_vector)
                        end
                else
                        println(this_vector)
                end

        end

end

all_mode_configuration(n,part.n_subset, only_photon_number_conserving  = true)
