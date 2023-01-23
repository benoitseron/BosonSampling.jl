include("packages_loop.jl")

function func(n)

    m = n

    n_lost = 1

    source = QuantumDot(13.5 / 80)

    d = Uniform(0,2pi)
    ϕ = nothing # rand(d,m)

    η_loss_lines = 0.9 * ones(m)
    η_loss_bs =  0.9 * ones(m-1)
    η_loss_source = get_η_loss_source(m,source)

    params = LoopSamplingParameters(n=n, η = η_thermalization(n), η_loss_bs = η_loss_bs, η_loss_lines = η_loss_lines, η_loss_source = η_loss_source, ϕ = ϕ)

    params_event = convert(SamplingParameters, params)
    params_event.o = ThresholdFockDetection(ThresholdModeOccupation(first_modes(n-n_lost, m).state))

    set_parameters!(params_event)

    ev = params_event.ev

    @time compute_probability!(ev)

end

a = func(20)

### previous method ###

params_no_η_loss_source = LoopSamplingParameters(n=n, η = η_thermalization(n), η_loss_bs = η_loss_bs, η_loss_lines = η_loss_lines, η_loss_source = nothing, ϕ = ϕ)

params_event = convert(SamplingParameters, params)
params_event.o = ThresholdFockDetection(ThresholdModeOccupation(first_modes(n-n_lost, m).state))

set_parameters!(params_event)

ev = params_event.ev

@time p_x_imperfect_source(params_event, 1., source)

### correcting the ThresholdDetection ###

n = 4
m = n

d = Uniform(0,2pi)
ϕ = nothing # rand(d,m)
η_loss_lines = nothing #0.2 * ones(m)
η_loss_bs = nothing # 0.2 * ones(m-1)

params = LoopSamplingParameters(n=n, η = η_thermalization(n), η_loss_bs = η_loss_bs, η_loss_lines = η_loss_lines, ϕ = ϕ)

params_event = convert(SamplingParameters, params)
params_event.o = ThresholdFockDetection(ThresholdModeOccupation([0,0,1,1]))

set_parameters!(params_event)

ev = params_event.ev
o = params_event.o

state = o.s.state

possible_threshold_detections(ev)

unique(possible_threshold_detections(n::Int, state::Vector{Int}))


# get all indexes with a photon

indexes = findall(x->x==1, state)

# get number of photons detected

n_detected = sum(state)

# finding the position of possible colliding photons
mode_configs_colliding_photons = all_mode_configurations(n - n_detected, n_detected, only_photon_number_conserving = false)

possible_states = []

# generate a list of each possible configuration of colliding photons
# add this configuration where a photon is detected in `state` (as labelled by `indexes`)

for mode_config in mode_configs_colliding_photons

    for i in 1:length(indexes)

        new_state = copy(state)

        new_state[indexes[i]] += mode_config[i]

        push!(possible_states, new_state)

    end

end

unique(possible_states)