using BosonSampling

isa_proba(0.1)

struct ThresholdDetectorDarkCount <: OutputMeasurementType
    s::ModeOccupation # no dark count output event
    dark_count_proba::Real
    ThresholdDetectorDarkCount(s::ModeOccupation, dark_count_proba::Real) = isa_proba(dark_count_proba) ? new() : error("$dark_count_proba is not a probability")

end

function compute_probability!(ev::Event{TIn,TOut}) where {TIn<:InputType, TOut<:ThresholdDetectorDarkCount}

    no_dark_count_event = Event(ev.input_state, FockDetection(ev.output_measurement.s), ev.interf)
    compute_probability!(no_dark_count_event)

    dark_count_proba_params = no_dark_count_event.proba_params

    deviation = zeros(ev.interferometer.m)

    ######## sample stochastically for each output mode
    # then transform s and output it

    # change to Bool.(s) for threshold detectors

end

is_collisionless([0,2])

wsample([0.1,0.9])
