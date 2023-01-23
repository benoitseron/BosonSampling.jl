# validation of experiments immediately from experimental data


"""

    function get_event_list(loaded_experiment::OneLoopData)

Takes in experimental data in the form of a `OneLoopData` and outputs a list of `Event` (as if detector readings) in a format ready to validate.

If loading from a file, do for instance 

    d = load("data/one_loop/n_6_m_12.jld")

    loaded_experiment = d["data"]

!!! Makes the following conversions by default 

    dict_sample_type_to_measurement = Dict(ThresholdModeOccupation => ThresholdFockDetection, ModeOccupation => FockDetection)

"""
function get_event_list(loaded_experiment::ExperimentalData)

    # we convert the experimental data into a list of event (as need for certifiers)

    @unpack samples, params, sample_type = loaded_experiment
    @unpack i, interferometer = params

    dict_sample_type_to_measurement = Dict(ThresholdModeOccupation => ThresholdFockDetection, ModeOccupation => FockDetection)

    if eltype(samples.proba) != Int
        error("this procedure requires an event count (registered as having an Int type)")
    end


    events = []
    for (count, number_observations) in zip(samples.counts, samples.proba)
        for j in 1:number_observations

            output_measurement = dict_sample_type_to_measurement[sample_type](count)
            event = Event(i, output_measurement, interferometer) 

            push!(events, event)
        end
    end

    events

end




