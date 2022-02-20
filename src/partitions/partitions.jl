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
        function Partition(subsets)
                check_subset_overlap(subsets)
                new(subsets, length(subsets))
        end
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
#
# struct PartitionCountingInterferometer <: Interferometer
#         partition_occupancy::PartitionOccupancy
#         physical_interferometer::Interferometer
#         fourier_indexes::Vector{Complex}
#         virtual_interferometer::Interferometer
#         # the phases are 2pi/(n+1) * fourier_indexes
#
#         function PartitionCountingInterferometer(partition_occupancy, physical_interferometer, fourier_indexes)
#
#                 """outputs the virtual interferometer that gives the fourier phase
#                 x(fourier_indexes), of which the multidimensional inverse
#                 fourier transform outputs the the probability to get partition_occupancy.counts
#                 photons in partition_occupancy.partition"""
#
#                 virtual_interferometer_matrix = copy(physical_interferometer.U)
#                 for i in length(partition_occupancy.partition.subsets)
#                         virtual_interferometer_matrix *= Diagonal(@. exp(2*pi*1im/(partition_occupancy.n+1) * fourier_indexes))
#                 end
#                 virtual_interferometer_matrix *= U'
#
#                 new(partition_occupancy, physical_interferometer, virtual_interferometer_matrix, fourier_indexes, virtual_interferometer_matrix)
#         end
#
# end

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


function compute_probabilities_partition(physical_interferometer::Interferometer, part::Partition, n::Int)

        """computes the probability to find a certain photon counts in a
        partition `part` of the output modes for the interferometer given

        returns : (counts = physical_indexes, probabilities = pdf)

        corresponding to the occupation numbers in the partition and the
        associated probability"""

        @warn "only implemented in the bosonic case so far"
        @warn "this gives what happens for input of type [1^n O^(m-n)]"

        @argcheck n <= m "more than one input per mode is not implemented"

        fourier_indexes = all_mode_configurations(n,part.n_subset, only_photon_number_conserving = false)
        probas_fourier = Array{ComplexF64}(undef, length(fourier_indexes))
        virtual_interferometer_matrix = similar(physical_interferometer.U)

        for (index_fourier_array, fourier_index) in enumerate(fourier_indexes)

                # for each fourier index, we recompute the virtual interferometer
                virtual_interferometer_matrix  = physical_interferometer.U
                diag = [one(eltype(virtual_interferometer_matrix)) for i in 1:m]

                for (i,fourier_element) in enumerate(fourier_index)
                        ##### this is where the wrong must be
                        ##### confusion vector dot operations
                        this_phase = exp(2*pi*1im/(n+1) * fourier_element)

                        @show fourier_element
                        @show this_phase
                        @show diag

                        for j in 1:length(diag)

                                if part.subsets[i].subset[j] == 1
                                        diag[j] *= this_phase

                                end
                                @show diag
                        end

                end

                virtual_interferometer_matrix *= Diagonal(diag)
                virtual_interferometer_matrix *= physical_interferometer.U'

                probas_fourier[index_fourier_array] = permanent(virtual_interferometer_matrix[1:n,1:n])
        end

        probas_fourier

        physical_indexes = copy(fourier_indexes)

        # probas_physical(physical_index) = 1/(n+1)^(part.n_subset) * sum([probas_fourier[i] * exp(-2*pi*1im/(n+1) * dot(physical_index, fourier_indexes[i])) for i in 1:length(fourier_indexes)])
        # previous function

        probas_physical(physical_index) = 1/(n+1)^(part.n_subset) * sum(probas_fourier[i] * exp(-2pi*1im/(n+1) * dot(physical_index, fourier_index)) for (i,fourier_index) in enumerate(fourier_indexes))


        pdf = [probas_physical(physical_index) for physical_index in physical_indexes]

        pdf = clean_pdf(pdf)

        (physical_indexes, pdf, probas_fourier)
end



function photon_number_conserving_events(physical_indexes, n; partition_spans_all_modes = false)

        """returns only the events conserving photon number n

        if partition_spans_all_modes = false, gives all events with less than n or n
        photons

        if partition_spans_all_modes = true only exact photon number conserving
        physical_indexes"""

        results = []
        for index in physical_indexes
                if partition_spans_all_modes == false
                        if sum(index) <= n
                                push!(results, index)
                        end
                else
                        if sum(index) == n
                                push!(results, index)
                        end
                end
        end
        results

end


function photon_number_non_conserving_events(physical_indexes,n ; partition_spans_all_modes = false)

        """returns the elements not conserving the number of photons"""

        setdiff(physical_indexes, photon_number_conserving_events(physical_indexes, n, ; partition_spans_all_modes = partition_spans_all_modes))

end

function check_photon_conservation(physical_indexes,  pdf, n; atol = ATOL, partition_spans_all_modes = false)

        """checks if probabilities corresponding to non photon number conserving
        events are zero"""

        events_to_check = photon_number_non_conserving_events(physical_indexes,n; partition_spans_all_modes = partition_spans_all_modes)

        for (i, index) in enumerate(physical_indexes)
                if index in events_to_check
                        @argcheck isapprox(clean_proba(pdf[i]),0, atol=atol)# "forbidden event has non zero probability"
                end
        end

end

function print_pdfs(physical_indexes, pdf, n; physical_events_only = false, partition_spans_all_modes = false)

        indexes_to_print = physical_events_only ? photon_number_conserving_events(physical_indexes, n; partition_spans_all_modes = partition_spans_all_modes) : physical_indexes

        println("---------------")
        println("Partition results : ")
        for (i, index) in enumerate(physical_indexes)
                if index in indexes_to_print
                        println("index = $index, p = $(pdf[i])")
                end

        end
        println("---------------")
end

### HOM tests: one mode ###

m = 2
n = 2

set1 = [1,0]
physical_interferometer = Fourier(m)
part = Partition([Subset(set1)])

(physical_indexes,  pdf) = compute_probabilities_partition(physical_interferometer, part, n)

photon_number_conserving_events(physical_indexes,n)

check_photon_conservation(physical_indexes, pdf, n)

print_pdfs(physical_indexes,  pdf, n)

### HOM tests: mode1, mode2 ###

m = 2
n = 2

set1 = [1,0]
set2 = [0,1]
physical_interferometer = Fourier(m)
part = Partition([Subset(set1), Subset(set2)])

(physical_indexes,  pdf) = compute_probabilities_partition(physical_interferometer, part, n)

print_pdfs(physical_indexes, pdf,n; partition_spans_all_modes = true, physical_events_only = true)

check_photon_conservation(physical_indexes, pdf, n; partition_spans_all_modes = true)






### other tests ###

m = 4
n = 3
set1 = zeros(Int,m)
set2 = zeros(Int,m)
set1[1:2] .= 1
set2[3:4] .= 1


physical_interferometer = RandHaar(m)
part = Partition([Subset(set1), Subset(set2)])

(physical_indexes,  pdf, probas_fourier) = compute_probabilities_partition(physical_interferometer, part, n)
fourier_indexes = copy(physical_indexes)

print_pdfs(physical_indexes,  pdf, n)
#print_pdfs(physical_indexes,  probas_fourier, n)

check_photon_conservation(physical_indexes, pdf, n; partition_spans_all_modes = false)
