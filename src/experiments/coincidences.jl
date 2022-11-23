"""

        get_fixed_number_coincidences(summarized_counts::MultipleCounts, number_coincidences::Int)

Takes in a `MultipleCount` containing event with different number of photons detected (when loss is present, that is) and returns a `MultipleCounts` with `number_coincidences` clicks in the detectors.

"""
function get_fixed_number_coincidences(summarized_counts::MultipleCounts, number_coincidences::Int)

        proba = Vector{Real}()
        counts = Vector{eltype(summarized_counts.counts)}()

        for i in 1:length(summarized_counts.proba)

                if typeof(summarized_counts.counts[1]) == ModeOccupation
                        state = summarized_counts.counts[i].state
                elseif typeof(summarized_counts.counts[1]) == ThresholdModeOccupation
                        state = summarized_counts.counts[i].clicks
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
