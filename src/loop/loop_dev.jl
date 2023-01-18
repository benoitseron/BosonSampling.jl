include("packages_loop.jl")

# making the formalism work for ThresholdModeOccupation
# basically this requires to compute the probabilities for each possible detected number of photons

n = 4
m = n

d = Uniform(0,2pi)
ϕ = nothing # rand(d,m)
η_loss_lines = nothing # 0.9 * ones(m)
η_loss_bs = nothing #1. * ones(m-1)

params = LoopSamplingParameters(n=n, η = η_thermalization(n), η_loss_bs = η_loss_bs, η_loss_lines = η_loss_lines, ϕ = ϕ)

params_event = convert(SamplingParameters, params)
params_event.o = ThresholdFockDetection(ThresholdModeOccupation([0,0,1,1]))

set_parameters!(params_event)

compute_probability!(params_event)

ev = params_event.ev

function compute_threshold_detection_probability(ev::Event{<:InputType,<:ThresholdFockDetection})

    state = ev.output_measurement.s.state

    detection_subset = Subset(state)
    complement = complement_subset(detection_subset)

    ##### need to check this stuff if loss

    part = BosonSampling.Partition([detection_subset, complement])

    # compute the probabilities associated with the partition created by the detection

    i = ev.input_state
    o = PartitionCountsAll(part)
    ev_partition = Event(i,o,ev.interferometer)

    @show mc = compute_probability!(ev_partition)

    # we are now interested in the case where there are zero photons in the complement, and from the number of threshold counts detected up to the number of input photons

    n_input = ev.input_state.r.n
    n_detected = ev.output_measurement.s.n_detected

    n_in_complement(i) = mc.counts[i].counts.state[2] 
    n_in_detectors(i) = mc.counts[i].counts.state[1]
    # there is potentially mc.counts[i].counts.state[3] lost

    function valid_count(i)
        if n_in_complement(i) == 0 && (n_in_detectors(i) >= n_detected && n_in_detectors(i) <= n_input)
            return true
        else
            return false
        end
    end

    probability_this_threshold_detection = 0

    for i in 1:length(mc.counts)
        if valid_count(i)
            probability_this_threshold_detection += mc.proba[i]
        end
    end

    probability_this_threshold_detection

end

compute_threshold_detection_probability(ev)


typeof(params)

n = 3
m = n

d = Uniform(0,2pi)
ϕ = nothing # rand(d,m)
η_loss_lines = 0.9 * ones(m)
η_loss_bs = 1. * ones(m-1)

params = LoopSamplingParameters(n=n, η = η_thermalization(n), η_loss_bs = η_loss_bs, η_loss_lines = η_loss_lines, ϕ = ϕ)

source = QuantumDot(efficiency = 0.85)

params_event = convert(SamplingParameters, params)
params_event_x = copy(params_event)

params_event_x.x = 0

params_event.o =  FockDetection(ModeOccupation([1,1,0]))
params_event_x.o =  FockDetection(ModeOccupation([1,1,0]))

set_parameters!(params_event)
set_parameters!(params_event_x)

compute_probability_imperfect_source(params_event, source)
compute_probability_imperfect_source(params_event_x, source)


p_x_imperfect_source(params_event, 0, source)

ev = params_event_x.ev
ev.output_measurement = FockDetection(ModeOccupation([1,1,1]))
ev

p_x_imperfect_source_update_this_event(ev, params_event_x, source)