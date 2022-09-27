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

function initialise_to_empty_vectors!(mc::MultipleCounts, type_proba, type_counts)

	mc.proba = Vector{type_proba}()
	mc.counts = Vector{type_counts}()

end

Base.show(io::IO, pb::MultipleCounts) = begin

	if pb.proba == nothing
		println(io, "Empty MultipleCounts")
	else
		for i in 1:length(pb.proba)

			println(io, "output: \n")
			println(io, pb.counts[i])
			println(io, "p = $(pb.proba[i])")
			println(io, "--------------------------------------")
		end
	end

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

Base.show(io::IO, ev::Event) = begin
	println("Event:\n")
	println("input state: ", ev.input_state.r, " (",get_parametric_type(ev.input_state)[1],")", "\n")
	println("output measurement: ", ev.output_measurement, "\n")
	println(ev.interferometer, "\n")
	println("proba_params: ", ev.proba_params)
end

struct GaussianEvent{TIn<:Gaussian, TOut<:OutputMeasurementType}

	input_state::GaussianInput{TIn}
	output_measurement::TOut
	interferometer::Union{Interferometer, Nothing}

	function GaussianEvent{TIn,TOut}(input_state, output_measurement, interferometer=nothing) where {TIn<:Gaussian, TOut<:OutputMeasurementType}
		new{TIn,TOut}(input_state, output_measurement, interferometer)
	end

	GaussianEvent(i,o,interf=nothing) = GaussianEvent{get_parametric_type(i)[1], get_parametric_type(o)[1]}(i,o,interf)

end

Base.show(io::IO, ev::GaussianEvent) = begin
	println("Event:\n")
	println("input state: ", ev.input_state.r, " (",get_parametric_type(ev.input_state)[1],")", "\n")
	println("output measurement: ", ev.output_measurement, "\n")
	println(ev.interferometer, "\n")
end


function check_probability_empty(ev::Event; resetting_message = true)
    if ev.proba_params.probability != nothing
		if resetting_message
			@warn "probability was already set in, rewriting"
		else
			@warn "unexpected probabilities found in Event"
		end
	end
end

Base.convert(::Type{Event{TIn, FockDetection}}, ev::Event{TIn, FockSample}) where {TIn <: InputType} = Event(ev.input_state, convert(FockDetection, ev.output_measurement), ev.interferometer, ev.proba_params)

# fs = FockSample([1,2,3])
# convert(FockDetection, fs)
