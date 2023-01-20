include("packages_loop.jl")

n = 8
m = n

n_lost = 4

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

params_no_η_loss_source = LoopSamplingParameters(n=n, η = η_thermalization(n), η_loss_bs = η_loss_bs, η_loss_lines = η_loss_lines, η_loss_source = nothing, ϕ = ϕ)

params_event = convert(SamplingParameters, params)
params_event.o = ThresholdFockDetection(ThresholdModeOccupation(first_modes(n-n_lost, m).state))

set_parameters!(params_event)

ev = params_event.ev

@time p_x_imperfect_source(params_event, 1., source)

