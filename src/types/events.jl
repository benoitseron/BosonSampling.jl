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
mutable struct Event{TIn<:InputType, TOut<:OutputMeasurementType}

    input_state::Input{TIn}
    output_measurement::TOut
    proba_params::EventProbability
    interferometer::Interferometer

    function Event{TIn,TOut}(input_state, output_measurement, interferometer::Interferometer, proba_params = nothing) where {TIn<:InputType, TOut<:OutputMeasurementType}

        if proba_params == nothing
            proba_params = EventProbability(nothing)
        end

		# in case input or outputs are given with m instead of 2m
		# for lossy cases, convert them to have the proper dimension for
		# computations

		if (LossParameters(typeof(interferometer)) == IsLossy())
			if input_state.m != 2*interferometer.m_real
				#println("converting Input to lossy")
				input_state = to_lossy(input_state)
			end

			if StateMeasurement(typeof(output_measurement)) == FockStateMeasurement()

				if output_measurement.s == nothing

					output_measurement.s = ModeOccupation(zeros(2*interferometer.m))

				elseif output_measurement.s.m != 2*interferometer.m_real
					#println("converting Output to lossy")
					to_lossy!(output_measurement)
				end
			end

		end

        new{TIn,TOut}(input_state, output_measurement, proba_params, interferometer)
    end

    Event(i,o,interf,p = nothing) = Event{get_parametric_type(i)[1], get_parametric_type(o)[1]}(i,o,interf,p)

end

# Base.show(io::IO, ev::Event) = begin
# 	println("Event:\n")
# 	println("input state: ", ev.input_state.r, " (",get_parametric_type(ev.input_state)[1],")", "\n")
# 	println("output measurement: ", ev.output_measurement, "\n")
# 	println(ev.interferometer, "\n")
# 	println("proba_params: ", ev.proba_params)
# end

struct GaussianEvent{TIn<:Gaussian, TOut<:OutputMeasurementType}

	input_state::GaussianInput{TIn}
	output_measurement::TOut
	interferometer::Union{Interferometer, Nothing}

	function GaussianEvent{TIn,TOut}(input_state, output_measurement, interferometer=nothing) where {TIn<:Gaussian, TOut<:OutputMeasurementType}
		new{TIn,TOut}(input_state, output_measurement, interferometer)
	end

	GaussianEvent(i,o,interf=nothing) = GaussianEvent{get_parametric_type(i)[1], get_parametric_type(o)[1]}(i,o,interf)

end

# Base.show(io::IO, ev::GaussianEvent) = begin
# 	println("Event:\n")
# 	println("input state: ", ev.input_state.r, " (",get_parametric_type(ev.input_state)[1],")", "\n")
# 	println("output measurement: ", ev.output_measurement, "\n")
# 	println(ev.interferometer, "\n")
# end


function check_probability_empty(ev::Event; resetting_message = true)
    if ev.proba_params.probability != nothing
		if resetting_message
			# @warn "probability was already set in, rewriting"
		else
			@warn "unexpected probabilities found in Event"
		end
	end
end

Base.convert(::Type{Event{TIn, FockDetection}}, ev::Event{TIn, FockSample}) where {TIn <: InputType} = Event(ev.input_state, convert(FockDetection, ev.output_measurement), ev.interferometer, ev.proba_params)

# fs = FockSample([1,2,3])
# convert(FockDetection, fs)

StateMeasurement(ev::Event) = StateMeasurement(typeof(ev.output_measurement))

# write possible_threshold_detections for an Event
# extract n from the input_state

function possible_threshold_detections(ev::Event)

    # check that the output_measurement is a ThresholdFockDetection

    @argcheck ev.output_measurement isa ThresholdFockDetection

    lossy = is_lossy(ev.interferometer)

    n = ev.input_state.r.n

    possible_threshold_detections(n,ev.output_measurement, lossy = lossy)

end

# write a function to deepcopy an Event importing the Base.deepcopy

function Base.copy(ev::Event)

	input_state = deepcopy(ev.input_state)
	output_measurement = deepcopy(ev.output_measurement)
	interferometer = deepcopy(ev.interferometer)
	proba_params = deepcopy(ev.proba_params)

	Event(input_state, output_measurement, interferometer, proba_params)

end

function to_threshold(ev::Event{TIn, FockDetection}) where {TIn <: InputType}

	input_state = deepcopy(ev.input_state)
	output_measurement = deepcopy(ev.output_measurement)
	interferometer = deepcopy(ev.interferometer)
	proba_params = deepcopy(ev.proba_params)

	output_measurement = convert(ThresholdFockDetection, output_measurement)

	Event(input_state, output_measurement, interferometer, proba_params)

end


n_clicks(ev::Event) = sum(ev.output_measurement.s.state)
max_possibly_lost_photons(ev::Event) = ev.input_state.n - n_clicks(ev)
