"""
	MultipleCounts()
	MultipleCounts(counts, proba)

Holds something like the photon counting probabilities with their respective
probability (in order to use them as a single observation). Can be declared
empty as a placeholder.

	Fields:
		- counts::Union{Nothing, Vector{ModeOccupation}, Vector{PartitionOccupancy}},
		- proba::Union{Nothing,Vector{Real}}
"""
mutable struct MultipleCounts

	counts::Union{Nothing, Vector{ModeOccupation}, Vector{PartitionOccupancy}}
	proba::Union{Nothing,Vector{Real}}

	MultipleCounts() = new(nothing,nothing)
	MultipleCounts(counts, proba) = new(counts,proba)

end

"""
	EventProbability(probability::Union{Nothing, Number})
	EventProbability(mc::MultipleCounts)

Holds the probability or probabilities of an [`Event`](@ref).

	Fields:
		- probability::Union{Number,Nothing, MultipleCounts}
		- precision::Union{Number,Nothing}
		- failure_probability::Union{Number,Nothing}
"""
mutable struct EventProbability
    probability::Union{Number,Nothing, MultipleCounts}
    precision::Union{Number,Nothing} # see remarks in conventions
    failure_probability::Union{Number,Nothing}

    function EventProbability(probability::Union{Nothing, Number})

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

	function EventProbability(mc::MultipleCounts)

		try
			mc.proba = clean_pdf(mc.proba)
			new(mc,nothing,nothing)
		catch
			error("invalid probability")
		end

	end

end

"""
	Event{TIn<:InputType, TOut<:OutputMeasurementType}

Event linking an input to an output.

	Fields:
		- input_state::Input{TIn}
		- output_measurement::TOut
		- proba_params::EventProbability
		- interferometer::Interferometer
"""
struct Event{TIn<:InputType, TOut<:OutputMeasurementType}

    input_state::Input{TIn}
    output_measurement::TOut
    proba_params::EventProbability
    interferometer::Interferometer

    function Event{TIn, TOut}(input_state, output_measurement, interferometer::Interferometer, proba_params = nothing) where {TIn<:InputType, TOut<:OutputMeasurementType}

        if proba_params == nothing
            proba_params = EventProbability(nothing)
        end

        new{TIn,TOut}(input_state, output_measurement, proba_params, interferometer)
    end

    Event(i,o,interf,p = nothing) = Event{get_parametric_type(i)[1], get_parametric_type(o)[1]}(i,o,interf,p)

end

function check_probability_empty(ev::Event)
    if ev.proba_params.probability != nothing
		@warn "probability was already set in, rewriting"
	end
end
