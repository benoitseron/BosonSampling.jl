function sample!(ev::Event{TIn, TOut}) where {TIn<:InputType, TOut <: FockSample}

    check_probability_empty(ev)

    if TIn == Distinguishable
        ev.output_measurement.s = ModeOccupation(classical_sampler(ev))
    elseif TIn == Bosonic
        ev.output_measurement.s = ModeOccupation(cliffords_sampler(ev))
    else
        error("not implemented")

    end

end



#
#
# n = 3
# m = 8
# interf = RandHaar(m)
# input_state = Input{Bosonic}(first_modes(n,m))
# o = FockSample()
#
# ev = Event(input_state, o, interf)
#
# sample!(ev)
