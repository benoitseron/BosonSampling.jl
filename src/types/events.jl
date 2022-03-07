mutable struct MultipleCounts

	"""holds something like the photon counting probabilities with their respective probability"""

	counts::Union{Nothing, Vector{ModeOccupation}, Vector{PartitionOccupancy}}
	proba::Union{Nothing,Vector{Real}}

	MultipleCounts() = new(nothing,nothing)

end

mutable struct EventProbability
    probability::Union{Number,Nothing, MultipleCounts}
    precision::Union{Number,Nothing} # see remarks in conventions
    failure_probability::Union{Number,Nothing}

    function EventProbability(probability = nothing)

        if probability == nothing
            new(nothing, nothing, nothing)
        else
            try
                probability = clean_proba(probability)
                new(probability,nothing,nothing)
            catch
                error("invalid probability")
            end
        end
    end
end


struct Event{TIn<:InputType, TOut<:OutputMeasurementType}

    input_state::Input{TIn}
    output_measurement::TOut
    proba_params::EventProbability
    interferometer::Interferometer

    function Event{TIn, TOut}(input_state, output_measurement, interferometer::Interferometer, proba_params = nothing) where {TIn<:InputType, TOut<:OutputMeasurementType}

        if proba_params == nothing
            proba_params = EventProbability()
        end

        new{TIn,TOut}(input_state, output_measurement, proba_params, interferometer)
    end

    Event(i,o,interf,p = nothing) = Event{get_parametric_type(i)[1], get_parametric_type(o)[1]}(i,o,interf,p)

end

function check_probability_empty(event::Event)
    if ev.proba_params.probability != nothing
		@warn "probability was already set in, rewriting"
	end
end
