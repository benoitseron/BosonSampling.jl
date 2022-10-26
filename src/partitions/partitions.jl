"""
    all_mode_configurations(n, n_subset; only_photon_number_conserving=false)
    all_mode_configurations(input_state::Input, part::Partition; only_photon_number_conserving=false)
    all_mode_configurations(input_state::Input, sub::Subset; only_photon_number_conserving=false)

Generate all possible photon counts of `n` photons in a partition/subset
of `n_subset` subsets.

!!! note
    - Does not take into account photon number conservation by default
    - This is the photon counting in partitions and not events outputs but it
      can be used likewise
"""
function all_mode_configurations(n,n_subset; only_photon_number_conserving = false)



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

all_mode_configurations(input_state::Input,part::Partition; only_photon_number_conserving = false) = all_mode_configurations(input_state.n,part.n_subset; only_photon_number_conserving = only_photon_number_conserving)

all_mode_configurations(input_state::Input,sub::Subset; only_photon_number_conserving = false) = all_mode_configurations(input_state.n,1; only_photon_number_conserving = only_photon_number_conserving)

"""
    remove_trivial_partitions!(part_list)

In a list of partitions sizes, ex. `[[2,0],[1,1],[0,2]]`, keeps only
the elements with non trivial subset size, in this ex. only `[1,1]`.
"""
function remove_trivial_partitions!(part_list)

    filter!(x -> !any(x .== 0), part_list)

end

"""
    ranked_partition_list(part_list)

Remove partitions such as `[1,2]` when `[2,1]` is already counted as only the
size of the partition counts; only keeps vectors of decreasing count.
"""
function ranked_partition_list(part_list)

    output_partitions = []

    for part in part_list

        keep_this_part = true
        for i in 1:length(part)-1
            if part[i] < part[i+1]
                keep_this_part = false
                break
            end
        end

        if keep_this_part
            push!(output_partitions, part)
        end
    end

    output_partitions
end

"""
    photon_number_conserving_events(physical_indexes, n; partition_spans_all_modes=false)

Return only the events conserving photon number `n`.

!!! note
    - If `partition_spans_all_modes`=`false`, gives all events with less than `n` or `n` photons
    - If `partition_spans_all_modes` = `true` only exact photon number conserving physical_indexes
"""
function photon_number_conserving_events(physical_indexes, n; partition_spans_all_modes = false)

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

"""
    photon_number_non_conserving_events(physical_indexes, n; partition_spans_all_modes=false)

Return the elements not conserving the number of photons.
"""
function photon_number_non_conserving_events(physical_indexes,n ; partition_spans_all_modes = false)

        setdiff(physical_indexes, photon_number_conserving_events(physical_indexes, n, ; partition_spans_all_modes = partition_spans_all_modes))

end

"""
    check_photon_conservation(physical_indexes,  pdf, n; atol=ATOL, partition_spans_all_modes=false)

Check if probabilities corresponding to non photon number conserving events are zero.
"""
function check_photon_conservation(physical_indexes,  pdf, n; atol = ATOL, partition_spans_all_modes = false)

        events_to_check = photon_number_non_conserving_events(physical_indexes,n; partition_spans_all_modes = partition_spans_all_modes)

        for (i, index) in enumerate(physical_indexes)
                if index in events_to_check
                        @argcheck isapprox(clean_proba(pdf[i]),0, atol=atol)# "forbidden event has non zero probability"
                end
        end

end

"""
    compute_probabilities_partition(physical_interferometer::Interferometer, part::Partition, input_state::Input)

Compute the probability to find a certain photon counts in a partition `part` of
the output modes for the given interferometer.

Return `(counts = physical_indexes, probabilities = pdf) corresponding to the
occupation numbers in the partition and the associated probability.
"""
function compute_probabilities_partition(physical_interferometer::Interferometer, part::Partition, input_state::Input)

        @argcheck at_most_one_photon_per_bin(input_state.r) "more than one input per mode is not implemented"

        n = input_state.n

        # if LossParameters(typeof(physical_interferometer)) == IsLossy()
        #         m = physical_interferometer.m_real
        # else
        m = input_state.m
        # end

        mode_occupation_list = fill_arrangement(input_state)
        S = input_state.G.S

        physical_indexes = []
        pdf = []

        if occupies_all_modes(part)
                # in this case we use a the trick of removing the last subset
                # and computing all as if in the partition without the last
                # subset and then to infer the photon number in the last one
                # form photon conservation

                (small_indexes, small_pdf) = compute_probabilities_partition(physical_interferometer, remove_last_subset(part), input_state)

                for (this_count, p) in zip(small_indexes, small_pdf)

                    new_count = copy(this_count)

                    if n-sum(this_count) >=0 # photon number conserving case

                            append!(new_count, n-sum(this_count))

                            push!(physical_indexes, new_count)
                            push!(pdf, p)
                    end

                end

        else

                fourier_indexes = all_mode_configurations(n,part.n_subset, only_photon_number_conserving = false)
                probas_fourier = Array{ComplexF64}(undef, length(fourier_indexes))
                virtual_interferometer_matrix = similar(physical_interferometer.U)

                for (index_fourier_array, fourier_index) in enumerate(fourier_indexes)

                        # for each fourier index, we recompute the virtual interferometer
                        virtual_interferometer_matrix  = physical_interferometer.U

                        diag = [1.0 + 0im for i in 1:m]
                        # this is not type stable
                        # but need it to be a complex float at least
                        for (i,fourier_element) in enumerate(fourier_index)

                                this_phase = exp(2*pi*1im/(n+1) * fourier_element)

                                for j in 1:length(diag)

                                        if part.subsets[i].subset[j] == 1

                                                diag[j] *= this_phase

                                        end

                                end

                        end

                        virtual_interferometer_matrix *= Diagonal(diag)
                        virtual_interferometer_matrix *= physical_interferometer.U'


                        # beware, only the modes corresponding to the
                        # virtual_interferometer_matrix[input_config,input_config]
                        # must be taken into account !
                        probas_fourier[index_fourier_array] = permanent(virtual_interferometer_matrix[mode_occupation_list,mode_occupation_list] .* S)
                end

                physical_indexes = copy(fourier_indexes)

                probas_physical(physical_index) = 1/(n+1)^(part.n_subset) * sum(probas_fourier[i] * exp(-2pi*1im/(n+1) * dot(physical_index, fourier_index)) for (i,fourier_index) in enumerate(fourier_indexes))


                pdf = [probas_physical(physical_index) for physical_index in physical_indexes]

        end

        pdf = clean_pdf(pdf)

        check_photon_conservation(physical_indexes, pdf, n; partition_spans_all_modes = occupies_all_modes(part))

        (physical_indexes, pdf)
end

"""
    compute_probability_partition_occupancy(physical_interferometer::Interferometer, part_occupancy::PartitionOccupancy, input_state::Input)

Compute the probability to find a partition occupancy.

!!! note
    Inefficient to use multiple times for the same physical setting, rather use
    compute_probabilities_partition.
"""
function compute_probability_partition_occupancy(physical_interferometer::Interferometer, part_occupancy::PartitionOccupancy, input_state::Input)

        (physical_indexes, pdf) = compute_probabilities_partition(physical_interferometer, part_occupancy.partition, input_state::Input)

        for (i,counts) in enumerate(physical_indexes)
                if counts == part_occupancy.counts.state
                        return pdf[i]
                end
        end
        nothing

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

"""
    compute_probability!(ev::Event{TIn,TOut}) where {TIn<:InputType, TOut<:PartitionCount}
    compute_probability!(ev::Event{TIn,TOut}) where {TIn<:InputType, TOut<:PartitionCountsAll}

Given a defined [`Event`](@ref), computes/updates its probability or set of probabilities
(for instance if looking at partition outputs, with `MultipleCounts` begin filled).

This function is defined separately as it is most often the most time consuming
step of calculations and one may which to separate the evaluation of probabilities
from preliminary definitions.
"""
function compute_probability!(ev::Event{TIn,TOut}) where {TIn<:InputType, TOut<:PartitionCount}

        check_probability_empty(ev)

        ev.proba_params.precision = eps()
        ev.proba_params.failure_probability = 0

        ev.proba_params.probability = compute_probability_partition_occupancy(ev.interferometer, ev.output_measurement.part_occupancy, ev.input_state)

end

function compute_probability!(ev::Event{TIn,TOut}) where {TIn<:InputType, TOut<:PartitionCountsAll}

        check_probability_empty(ev)



        ev.proba_params.precision = eps()
        ev.proba_params.failure_probability = 0

        i = ev.input_state
        part = ev.output_measurement.part

        if i.m != part.m
                if i.m == 2*part.m
                        # @warn "converting the partition to a lossy one"
                        part = to_lossy(part)
                        ev.output_measurement = PartitionCountsAll(part)
                else
                        error("incompatible i, part")
                end
        end

        (part_occ,  pdf) = compute_probabilities_partition(ev.interferometer, ev.output_measurement.part, i)

        # clean up to keep photon number conserving events (possibly lossy events in the partition occupies all modes)

        part_occ_physical = []
        pdf_physical = []

        n = i.n

        for (i,occ) in enumerate(part_occ)
            if sum(occ) <= n
                    if occupies_all_modes(ev.output_measurement.part)
                            if sum(occ) == n

                                push!(part_occ_physical, occ)
                                push!(pdf_physical, pdf[i])

                            end
                    else
                            push!(part_occ_physical, occ)
                            push!(pdf_physical, pdf[i])
                    end

            end
        end

        mc = MultipleCounts([PartitionOccupancy(ModeOccupation(occ),n,part) for occ in part_occ_physical], pdf_physical)

        ev.proba_params.probability = EventProbability(mc).probability

end


"""
    to_partition_count(event::Event{TIn, TOut}, part::Partition) where {TIn<:InputType, TOut <: FockDetection}

Converts an `Event` with `FockDetection` to a `PartitionCount` one.
"""
function to_partition_count(ev::Event{TIn, TOut}, part::Partition) where {TIn<:InputType, TOut <: Union{FockDetection, FockSample}}

    n_subsets = part.n_subset

    counts_array = zeros(Int,n_subsets)

    for i in 1:n_subsets
        counts_array[i] = sum(ev.output_measurement.s.state .* part.subsets[i].subset)
    end

    @argcheck sum(counts_array) == ev.input_state.n

    o = PartitionCount(PartitionOccupancy(ModeOccupation(counts_array), ev.input_state.n, part))

    new_ev = Event(ev.input_state, o, ev.interferometer)

    new_ev

end

"""
    p_partition(ev::Event{TIn1, TOut1}, ev_theory::Event{TIn2, TOut2}) where {TIn1<:InputType, TOut1 <: PartitionCount, TIn2 <:InputType, TOut2 <:PartitionCountsAll}

Outputs the probability that an observed count in `ev` happens under the conditions set by `ev_theory`.
For instance, if we take the conditions

    ib = Input{Bosonic}(first_modes(n,m))
    part = equilibrated_partition(m,n_subsets)
    o = PartitionCountsAll(part)

    evb = Event(ib,o,interf)

then

    p_partition(ev, evb)

gives the probability that this `ev` is observed under the hypotheses of `ev_theory`.
"""
function p_partition(ev::Event{TIn1, TOut1}, ev_theory::Event{TIn2, TOut2}) where {TIn1<:InputType, TOut1 <: PartitionCount, TIn2 <:InputType, TOut2 <:PartitionCountsAll}

        #check interferometer, input configuration
    @argcheck ev.output_measurement.part_occupancy.partition == ev_theory.output_measurement.part
    @argcheck ev.interferometer == ev_theory.interferometer
    @argcheck ev.input_state.r == ev_theory.input_state.r

    # compute the probabilities if they are not already known
    ev_theory.proba_params.probability == nothing ? compute_probability!(ev_theory) : nothing

    p = ev_theory.proba_params.probability

    observed_count = ev.output_measurement.part_occupancy.counts.state

    # look up the probability of this count

    proba_this_count = nothing

    for (proba, theoretical_count) in zip(p.proba, p.counts)
        if observed_count == theoretical_count.counts.state
            proba_this_count = proba
            break
        end
    end

    proba_this_count

end

function compute_probability!(params::PartitionSamplingParameters)

    @unpack n, m, interf, T, mode_occ, i, n_subsets, part, o, ev = params

    compute_probability!(ev)

end

"""
partition_thermalization(m)

Defines the last mode, single mode subset for thermalization. This corresponds to the first mode of the interferometer with spatial bins (to be checked).
"""
partition_thermalization(m) = begin

s1 = Subset(ModeList(m,m))
s2 = Subset(first_modes(m-1,m))
Partition([s1,s2])

end
