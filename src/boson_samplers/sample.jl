function sample!(ev::Event{TIn, TOut}) where {TIn<:InputType, TOut <: FockSample}

    check_probability_empty(ev)

    if TIn == Distinguishable
        ev.output



n = 3
m = 8
interf = RandHaar(m)
input_state = Input{Bosonic}(first_modes(n,m))
o

ev = Event(input_state, FockSample())
