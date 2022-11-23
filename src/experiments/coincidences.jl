"""

        get_fixed_number_coincidences(summarized_counts::MultipleCounts, number_coincidences::Int)

Takes in a `MultipleCount` containing event with different number of photons detected (when loss is present, that is) and returns a `MultipleCounts` with `number_coincidences` clicks in the detectors.

"""
function get_fixed_number_coincidences(summarized_counts::MultipleCounts, number_coincidences::Int)

        proba = Vector{Real}()
        counts = Vector{ModeOccupation}()

        for i in 1:length(summarized_counts.proba)

                if sum(summarized_counts.counts[i].state) == number_coincidences
                        push!(proba, summarized_counts.proba[i])
                        push!(counts, summarized_counts.counts[i])
                end
        end

        MultipleCounts(counts, proba)
end
