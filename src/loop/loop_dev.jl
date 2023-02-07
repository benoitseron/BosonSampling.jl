include("packages_loop.jl")

### verifying threshold probabilities ###

n = 4
m = 2n

interf = Fourier(m)

state = [0,0,1,1]
o = ThresholdFockDetection(ThresholdModeOccupation(state))
i = Input{Bosonic}(first_modes(n,m))

ev = Event(i, o, interf)

compute_probability!(ev)


### full threshold distribution ###

o = BosonSamplingThresholdDistribution()
ev_full_distribution_threshold = Event(i, o, interf)

compute_probability!(ev_full_distribution_threshold)

o = BosonSamplingDistribution()
ev_full_distribution = Event(i, o, interf)

compute_probability!(ev_full_distribution)

mc = ev_full_distribution.proba_params.probability
mc_thresholdised = to_threshold(mc)

mc_threshold = ev_full_distribution_threshold.proba_params.probability

sort!(mc_thresholdised) 
sort!(mc_threshold)

mc_thresholdised.proba ≈ mc_threshold.proba

deepcopy(mc_thresholdised)






ev_full_distribution.proba_params.probability

ev_full_distribution
# get all indexes with a photon

indexes = findall(x->x==1, state)

# get number of photons detected

n_detected = sum(state)

# if no loss, nothing to do 
if n_detected == n
    return [state]
end

# finding the position of possible colliding photons
mode_configs_colliding_photons = all_mode_configurations(n - n_detected, n_detected, only_photon_number_conserving = true)

possible_states = []

# generate a list of each possible configuration of colliding photons
# add this configuration where a photon is detected in `state` (as labelled by `indexes`)



for mode_config in mode_configs_colliding_photons

    # @show mode_config
 
    new_state = copy(state)

        # @show new_state
        # @show indexes

    for i in 1:length(indexes)

        new_state[indexes[i]] += mode_config[i]

    end

        # @show new_state

        push!(possible_states, new_state)

end

unique(possible_states)

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

ev_

# Event{Bosonic, FockDetection}(Input{Bosonic}(state = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 10, 20, GramMatrix{Bosonic}(10, ComplexF64[1.0 + 0.0im 1.0 + 0.0im … 1.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 1.0 + 0.0im … 1.0 + 0.0im 1.0 + 0.0im; … ; 1.0 + 0.0im 1.0 + 0.0im … 1.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 1.0 + 0.0im … 1.0 + 0.0im 1.0 + 0.0im], nothing, nothing, OrthonormalBasis(nothing)), nothing), FockDetection(state = [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0]), EventProbability(2.4682721085522713e-48, 2.220446049250313e-16, 0), Interferometer :

# Type : LossyCircuit
# m : 20
# U : 
# Complex[0.12656250000000002 + 0.0im -0.06200979635307635 + 0.0im … 0.0 + 0.0im 0.0 + 0.0im; 0.11161763343553743 + 0.0im 0.07031250000000001 + 0.0im … 0.0 + 0.0im 0.0 + 0.0im; … ; 0.0 + 0.0im 0.0 + 0.0im … 0.16875 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im … 0.0 + 0.0im 0.16875 + 0.0im])



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

# compute_probability!(ev)

ev.input_state.r

i = Input{Bosonic}(ev.input_state.r)

possible_threshold_detections(ev)

state = ev.output_measurement.s.state
n = ev.input_state.r.n

possible_threshold_detections(n, state, lossy = true)

possible_outputs = possible_threshold_detections(ev)

# REMOVING THE PARTIAL DISTINGUISHABILITY

i = Input{Bosonic}(ev.input_state.r)
ev_ = Event(i, FockDetection(possible_outputs[1]), ev.interferometer)



# for each possible output, compute the probability
# and sum them up

# Event{OneParameterInterpolation, FockDetection}(Input{OneParameterInterpolation}(state = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 10, 20, GramMatrix{OneParameterInterpolation}(10, [1.0 0.9 … 0.9 0.9; 0.9 1.0 … 0.9 0.9; … ; 0.9 0.9 … 1.0 0.9; 0.9 0.9 … 0.9 1.0], nothing, 0.9, OrthonormalBasis(nothing)), 0.9), FockDetection(state = [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0]), EventProbability(nothing, nothing, nothing), Interferometer :

# Type : LossyCircuit
# m : 20
# U : 
# Complex[-0.13636899171668213 + 0.0im 0.10406410456750137 + 0.0im … 7.092747648714338e-5 + 0.0im 0.0 + 0.0im; -0.13636899171668215 + 0.0im -0.38901677063462753 + 0.0im … -0.00026514404719058915 + 0.0im 0.0 + 0.0im; … ; 0.0 + 0.0im 0.0 + 0.0im … -0.5092793132753864 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im … 0.0 + 0.0im -0.19285487757287617 + 0.0im])


@time compute_probability!(ev_)

p = 0.
@showprogress for o in possible_outputs
    ev_.output_measurement = FockDetection(o)
    p += compute_probability!(ev_)
end

ev_ = Event(ev.input_state, FockDetection(possible_outputs[1]), ev.interferometer)