include("packages_loop.jl")

### correcting the ThresholdDetection ###

n = 10
m = n

d = Uniform(0,2pi)
ϕ = nothing # rand(d,m)
η_loss_lines = 0.9 * ones(m)
η_loss_bs = nothing # 0.2 * ones(m-1)
η_loss_source = get_η_loss_source(m,QuantumDot(13.5 / 80))

params = LoopSamplingParameters(n=n, η = η_thermalization(n), η_loss_bs = η_loss_bs, η_loss_lines = η_loss_lines, η_loss_source = η_loss_source, ϕ = ϕ)

params_event = convert(SamplingParameters, params)
params_event.o = ThresholdFockDetection(ThresholdModeOccupation([0,0,0,1,1,1,1,1,1,1]))

set_parameters!(params_event)

ev = params_event.ev

state = ev.output_measurement.s.state
n = ev.input_state.r.n

possible_outputs = possible_threshold_detections(ev)

ev_ = Event(ev.input_state, FockDetection(possible_outputs[1]), ev.interferometer)

@time compute_probability!(ev_)









compute_probability!(params_event)

### too large permanents computed ###

d = load("/home/benoitseron/Documents/thesis/one_loop_sampling/code/data/one_loop/n_10_m_10.jld")

loaded_experiment = d["data"]

@unpack params, samples = loaded_experiment

(η_lines, η_bs) = (0.8624409073596869, 0.9346111367443221)
# estimated parameters

set_uniform_losses!(η_lines, η_bs, params)

samples

events = get_event_list(loaded_experiment)

events[1].input_state

# ##### REDUCING EVENT LIST, to be removed #####

events = events[1:10]

ev = events[1]

compute_probability!(ev)

ev.input_state 

possible_threshold_detections(ev)

state = ev.output_measurement.s.state
n = ev.input_state.r.n

possible_threshold_detections(n, state, lossy = true)

possible_outputs = possible_threshold_detections(ev)

ev_ = Event(ev.input_state, FockDetection(possible_outputs[1]), ev.interferometer)

# for each possible output, compute the probability
# and sum them up

@time compute_probability!(ev_)

p = 0.
@showprogress for o in possible_outputs
    ev_.output_measurement = FockDetection(o)
    p += compute_probability!(ev_)
end

ev_ = Event(ev.input_state, FockDetection(possible_outputs[1]), ev.interferometer)