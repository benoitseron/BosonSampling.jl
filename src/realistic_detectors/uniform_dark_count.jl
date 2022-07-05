using BosonSampling

struct ThresholdDetectorDarkCount <: OutputMeasurementType
    s_no_dark_count::ModeOccupation # no dark count output event
    s::Union{ModeOccupation, Nothing} # observed event, set as nothing by default
    dark_count_proba::Real

    ThresholdDetectorDarkCount(s_no_dark_count::ModeOccupation, dark_count_proba::Real) = isa_proba(dark_count_proba) ? new(s_no_dark_count, nothing, dark_count_proba) : error("$dark_count_proba is not a probability")

end

function compute_probability!(ev::Event{TIn,TOut}) where {TIn<:InputType, TOut<:ThresholdDetectorDarkCount}

    no_dark_count_event = Event(ev.input_state, FockDetection(ev.output_measurement.s_no_dark_count), ev.interf)
    compute_probability!(no_dark_count_event)

    dark_count_proba_params = no_dark_count_event.proba_params

    sample_dark_count() = wsample(p,1-p) == 1 ? 1 : 0
    deviation = ModeOccupation([sample_dark_count() for i in 1:ev.interferometer.m])

    output_state = ev.output_measurement.s_no_dark_count + deviation
    # dark counts added on top of the normal counts

    # change to boolean values for threshold detectors
    ev.output_measurement.s = ModeOccupation(@. !iszero(state))
end

#### the thing is that it is not really the actual probability of this event... we should montecarlo and inverse to get this specific event
#### or do something more akin to a sampling function

n = 2
m = 4
i = Input{Bosonic}(first_modes(n,m))
o = FockDetection([1,0,1,0])
interf = SingleModeNonLinearPhaseShift(m, π)
