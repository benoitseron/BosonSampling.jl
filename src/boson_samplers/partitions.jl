# function sample!(ev::Event{TIn, TOut}) where {TIn<:InputType, TOut <: PartitionSample}
#
#     check_probability_empty(ev)
#
#     if TIn == Distinguishable
#         ev.output_measurement.s = ModeOccupation(classical_sampler(ev))
#     elseif TIn == Bosonic
#         ev.output_measurement.s = ModeOccupation(cliffords_sampler(ev))
#     else
#         error("not implemented")
#     end
#
# end
#
#
# ev = Event(ib,PartitionCount(wsample(pb.counts, pb.proba)), interf)


########this would not be efficient because recomputing compute_probability!(ev_theory) each time...
####### although we could store it somewhere and pass it out
