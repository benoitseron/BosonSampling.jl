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
