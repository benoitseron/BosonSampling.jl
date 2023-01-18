include("packages_loop.jl")



vtypeof(params)

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