"""
    sample!(ev::Event{TIn, TOut}) where {TIn<:InputType, TOut <: FockSample}
    sample!(ev::Event{TIn,TOut}, loss::Real) where {TIn<:InputType, TOut<:FockSample}

Simulate a boson sampling experiment form a given [`Event`](@ref).
"""
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

function sample!(ev::Event{TIn,TOut}, loss::Real) where {TIn<:InputType, TOut<:FockSample}

    check_probability_empty(ev)

    if TIn <: PartDist
        if TIn == OneParameterInterpolation
            ev.output_measurement.s = ModeOccupation(noisy_sampler(ev,loss))
        else
            error("model is not valid")
        end
    else
        error("not implemented")
    end
end

"""
    scattershot_sampling(n::Int, m::Int; N=1000, interf=nothing)

Simulate `N` times a scattershot boson sampling experiment with `n` photons among `m` modes.
The interferometer is set to [`RandHaar`](@ref) by default.
"""
function scattershot_sampling(n::Int, m::Int; N=1000, interf=nothing)

    out_ = Vector{Vector{Int64}}(undef, N)
    in_ = Vector{Vector{Int64}}(undef, N)

    for i in 1:N
        input = Input{Bosonic}(ModeOccupation(random_occupancy(n,m)))
        interf == nothing ? F = RandHaar(m) : F = interf
        ev = Event(input, FockSample(), F)
        sample!(ev)

        in_[i] = input.r.state
        out_[i] = ev.output_measurement.s.state
    end

    joint = sort([[in_[i],out_[i]] for i in 1:N])
    res = counter(joint)

    k = collect(keys(res))
    k = reduce(hcat,k)'
    vals = collect(values(res))

    k1 = unique(k[:,1])
    k2 = unique(k[:,2])
    M = zeros(length(k1),length(k2))
    for i in 1:length(k1)
        for j in 1:length(k2)
            M[i,j] = res[[k1[i]',k2[j]']]
        end
    end

    k1 = [string(i) for i in k1]
    k2 = [string(i) for i in k2]
    heatmap(k1,k2,M')
end
