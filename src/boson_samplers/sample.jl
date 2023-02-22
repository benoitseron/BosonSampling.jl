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
        ev.output_measurement.s = ModeOccupation(clifford_sampler_unoptimised(ev))
    elseif TIn == OneParameterInterpolation
        ev.output_measurement.s = ModeOccupation(noisy_sampler(ev,1))
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

# define the sampling algorithm:
# the function sample! is modified to take into account
# the new measurement type
# this allows to keep the same syntax and in fact reuse
# any function that would have previously used sample!
# at no cost
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

function sample!(ev::Event{TIn, TOut}) where {TIn<:InputType, TOut <: RealisticDetectorsFockSample}

    # sample with dark counts but seeing all photons
    ev_dark = Event(ev.input_state, DarkCountFockSample(ev.output_measurement.p_dark), ev.interferometer)
    sample!(ev_dark)
    sample_dark = ev_dark.output_measurement.s

    # remove each of the readings with p_no_count
    for mode in 1:sample_dark.m
        if do_with_probability(ev.output_measurement.p_no_count)
            sample_dark.state[mode] = 0
        end
    end

    ev.output_measurement.s = sample_dark
end

function sample!(params::SamplingParameters)
    sample!(params.ev)
end

function get_lossy_threshold_detector_reading(sampling_params::SamplingParameters)

        this_sample = BosonSampling.sample!(sampling_params)
        this_sample = ModeOccupation(this_sample.state[1:sampling_params.n]) #only keeping the relevant modes, not the environment

        sum(this_sample.state) # number of remaining photons

        to_threshold(this_sample)

end

function generate_approximate_threshold_lossy_distribution(sampling_params::SamplingParameters, n_samples::Int)

        samples = Vector{ModeOccupation}()
        counts = Vector{Int}()

        @showprogress for sample_count in 1:n_samples

                this_sample = get_lossy_threshold_detector_reading(sampling_params)

                index_array = findall(x->x==this_sample, samples)

                if length(index_array) == 0
                        push!(samples, this_sample)
                        push!(counts, 1)
                elseif length(index_array) == 1
                        counts[index_array[1]] += 1
                end
        end

        MultipleCounts(samples, counts ./ sum(counts)) # normalize to make it a probability as is eaten by MultipleCounts

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
