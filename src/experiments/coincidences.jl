"""

        get_fixed_number_coincidences(summarized_counts::MultipleCounts, number_coincidences::Int)

Takes in a `MultipleCount` containing event with different number of photons detected (when loss is present, that is) and returns a `MultipleCounts` with `number_coincidences` state in the detectors.

"""
function get_fixed_number_coincidences(summarized_counts::MultipleCounts, number_coincidences::Int)

        proba = Vector{Real}()
        counts = Vector{eltype(summarized_counts.counts)}()

        for i in 1:length(summarized_counts.proba)

                if typeof(summarized_counts.counts[1]) == ModeOccupation
                        state = summarized_counts.counts[i].state
                elseif typeof(summarized_counts.counts[1]) == ThresholdModeOccupation
                        state = summarized_counts.counts[i].state
                end

                if sum(state) == number_coincidences
                        push!(proba, summarized_counts.proba[i])
                        push!(counts, summarized_counts.counts[i])
                end
        end

        if length(proba) == 0
                @warn "empty counts for number_coincidences = $number_coincidences"
        end

        MultipleCounts(counts, proba)
end


"""
    max_coincidences(samples::MultipleCounts)

Gives the maximum number of clicks observed.
"""
max_coincidences(samples::MultipleCounts) = maximum([sum(samples.counts[i].state) for i in 1:length(samples.counts)])


"""

    sort_by_coincidence_counts(samples::MultipleCounts)

Outputs a dictionary linking the number of coincidences observed to the relevant samples.

"""
function sort_by_coincidence_counts(samples::MultipleCounts)

    sorted_samples = Dict{Int, MultipleCounts}()

    n_max = max_coincidences(samples)

    for n_coincidences in 1:n_max

        sorted_samples[n_coincidences] = get_fixed_number_coincidences(samples, n_coincidences)

    end

    sorted_samples

end


function get_fixed_number_photons(mc::MultipleCounts, n_photons_detected::Int)

        proba = Vector{Real}()
        counts = Vector{eltype(mc.counts)}()

        for (p, count) in zip(mc.proba, mc.counts)

                if sum(count) == n_photons_detected
                        push!(proba, p)
                        push!(counts, count)
                end

        end

        MultipleCounts(counts, proba)

end

function sort_by_detected_photons(samples::MultipleCounts)

        sorted_samples = Dict{Int, MultipleCounts}()
    
        n_max = max_coincidences(samples)
    
        for n_coincidences in 1:n_max
    
            sorted_samples[n_coincidences] = get_fixed_number_photons(samples, n_coincidences)
    
        end
    
        sorted_samples
    
end