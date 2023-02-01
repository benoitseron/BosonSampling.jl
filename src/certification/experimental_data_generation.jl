function generate_experimental_data(;n_events, n, m, interf, TIn, x = nothing, input_mode_occupation = nothing)

    if input_mode_occupation == nothing
        input_mode_occupation = first_modes(n,m)
    end

    if TIn == OneParameterInterpolation

        @argcheck x != nothing "need to set partial distinguishability"

        input_state = Input{TIn}(input_mode_occupation, x)

    elseif TIn == Bosonic || TIn == Distinguishable

        input_state = Input{TIn}(input_mode_occupation)
    end

    events = []

    for i in 1:n_events

        # note that we don't compute the event probability
        # as we would just have experimental observations
        # of counts

        ev = Event(input_state, FockSample(), interf)
        BosonSampling.sample!(ev)

        ev = convert(Event{TIn, FockDetection}, ev)

        push!(events, ev)

    end

    if !isa(events, Vector{Event})
        events = convert(Vector{Event}, events)
    end

    events
end

function generate_experimental_data(params::LoopSamplingParameters, n_events::Int)

    @unpack n, m, interferometer, T, x, i = params


    generate_experimental_data(;n_events = n_events, n = n, m = m, interf = interferometer, TIn = T, x = x, input_mode_occupation = i.r)
end

function generate_experimental_data(loaded_experiment::OneLoopData, n_events::Int)

    @unpack n, m, interferometer, T, x, i = loaded_experiment.params

    data = generate_experimental_data(loaded_experiment.params, n_events)

    if loaded_experiment.sample_type == ThresholdModeOccupation

        data_threshold = []
        println("converting to threshold mode occupation")
        for ev in data
            push!(data_threshold, to_threshold(ev))
            
        end
        return data_threshold
    end

    data
end