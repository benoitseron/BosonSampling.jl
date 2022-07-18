"""
    classical_sampler(U, n, m)
    classical_sampler(;input::Input, interf::Interferometer)

Sample photons according to the [`Distinguishable`](@ref) case.
"""
function classical_sampler(U, n, m)

    output_state = zeros(Int,m)
    output_modes = collect(1:m)

    #@warn "check U or U'"
    for j in 1:n
        this_output_mode = wsample(output_modes, abs.(U[j,:]) .^2)
        output_state[this_output_mode] += 1
    end

    output_state

end

classical_sampler(;input::Input, interf::Interferometer) = classical_sampler(interf.U, input.n, input.m)

"""
    classical_sampler(ev::Event{TIn, TOut}; occupancy_vector = true) where {TIn<:InputType, TOut <: FockSample}

Sampler for an `Event`. Note the difference of behaviour if `occupancy_vector = true`.
"""
function classical_sampler(ev::Event{TIn, TOut}; occupancy_vector = true) where {TIn<:InputType, TOut <: FockSample}
    s = classical_sampler(input = ev.input_state, interf = ev.interferometer)

    occupancy_vector ? mode_occupancy_to_occupancy_vector(s, ev.input_state.m) : s
end
