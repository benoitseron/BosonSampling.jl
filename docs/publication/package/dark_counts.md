```julia
# define a new measurement type of a simple dark counting detector

mutable struct DarkCountFockSample <: OutputMeasurementType

    s::Union{ModeOccupation, Nothing} # observed output, possibly undefined
    p::Real # probability of a dark count in each mode

    DarkCountFockSample(p::Real) = isa_probability(p) ? new(nothing, p) : error("invalid probability")
    # instantiate if no known output
end

# works DarkCountFockSample(0.01)
# fails DarkCountFockSample(-1)

# define the sampling algorithm
function sample!(ev::Event{TIn, TOut}) where {TIn<:InputType, TOut <: DarkCountFockSample}


    # sample without dark counts
    ev_no_dark = Event(ev.input_state, FockSample(), ev.interferometer)
    sample!(ev_no_dark)
    sample_no_dark = ev_no_dark.output_measurement.s

    # now, apply the dark counts to "perfect" samples

    observe_dark_count(p) = Int(do_with_probability(p)) # 1 with probability p, 0 with probability 1-p
    dark_counts = [observe_dark_count(ev.output_measurement.p) for i in 1: ev.input_state.m]

    ev.output_measurement.s = sample_no_dark + dark_counts
end

n = 10
m = 10
p_dark = 0.1
# experiment parameters
input_state = first_modes(n,m)
interf = RandHaar(m)
i = Input{Bosonic}(input_state)
o = DarkCountFockSample(p_dark)
ev = Event(i,o,interf)

sample!(ev)

# DOESNT WORK but ok if I execute the code of sample copy pasted below
# n = 3
# m = 8
# interf = RandHaar(m)
# input_state = Input{Bosonic}(first_modes(n,m))
# o = FockSample()
#
# ev = Event(input_state, o, interf)
#
# sample!(ev)

#
# function sample!(ev::Event{TIn, TOut}) where {TIn<:InputType, TOut <: FockSample}
#
#     check_probability_empty(ev)
#
#     if TIn == Distinguishable
#         ev.output_measurement.s = ModeOccupation(classical_sampler(ev))
#     elseif TIn == Bosonic
#         ev.output_measurement.s = ModeOccupation(cliffords_sampler(ev))
#     else
#         error("not implemented")
#
#     end
#
# end
#
# function do_with_probability(p)
#
#     """returns true with probability p, false with (1-p)"""
#
#     rand() < p ? true : false
#
# end
```