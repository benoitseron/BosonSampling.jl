# define a new measurement type of a simple dark counting detector

mutable struct DarkCountFockSample <: OutputMeasurementType

    s::Union{ModeOccupation, Nothing} # observed output, possibly undefined
    p::Real # probability of a dark count in each mode

    DarkCountFockSample(p::Real) = isa_probability(p) ? new(nothing, p) : error("invalid probability")
    # instantiate if no known output
end

# works DarkCountFockSample(0.01)
# fails DarkCountFockSample(-1)

# define the sampling algorithm:
# the function sample! is modified to take into account
# the new measurement type
# this allows to keep the same syntax and in fact reuse
# any function that would have previously used sample!
# at no cost
function BosonSampling.sample!(ev::Event{TIn, TOut}) where {TIn<:InputType, TOut <: DarkCountFockSample}

    # sample without dark counts
    ev_no_dark = Event(ev.input_state, FockSample(), ev.interferometer)
    sample!(ev_no_dark)
    sample_no_dark = ev_no_dark.output_measurement.s

    # now, apply the dark counts to "perfect" samples
    observe_dark_count(p) = Int(do_with_probability(p)) # 1 with probability p, 0 with probability 1-p
    dark_counts = [observe_dark_count(ev.output_measurement.p) for i in 1: ev.input_state.m]

    ev.output_measurement.s = sample_no_dark + dark_counts
end

### example ###

# experiment parameters
n = 10
m = 10
p_dark = 0.1
input_state = first_modes(n,m)
interf = RandHaar(m)
i = Input{Bosonic}(input_state)
o = DarkCountFockSample(p_dark)
ev = Event(i,o,interf)

sample!(ev)
# output:
# state = [3, 1, 0, 3, 0, 1, 2, 0, 0, 0]
