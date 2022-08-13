# define a new measurement type of a simple dark counting detector

struct DarkCountFockSample <: OutputMeasurementType

    s::Union{ModeOccupation, Nothing} # observed output, possibly undefined
    p::Real # probability of a dark count in each mode

    DarkCountFockSample(p::Real) = isa_probability(p) ? new(nothing, p) : error("invalid probability")
    # instantiate if no known output
end

# works DarkCountFockSample(0.01)
# fails DarkCountFockSample(-1)

# convert dark-count events to no-dark-count events for sampling

function convert(Event{TIn, FockSample}, ev::Event{TIn, DarkCountFockSample}) where {TIn <: InputType}

    Event{TIn, FockSample}(ev.input_state, FockSample(), ev.proba_params, ev.interferometer)

end


# define the sampling algorithm
function sample!(ev::Event{TIn, TOut}) where {TIn<:InputType, TOut <: DarkCountFockSample}


    # sample without dark counts
    ev_no_dark = convert(Event{TIn, FockSample}, ev)
    sample!(ev_no_dark)
    sample_no_dark = ev_no_dark.output_measurement.s

    # now, apply the dark counts to "perfect" samples

    observe_dark_count(p) = Int(do_with_probability(p)) # 1 with probability p, 0 with probability 1-p
    dark_counts = [observe_dark_count(ev.output_measurement.p) for i in 1: ################# lacking a definition of the photon number!]

end

isa_probability(0.5)

isa_probability(-1)
