include("packages_loop.jl")

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



"""
    p_x_imperfect_source(params_event::SamplingParameters ,x, source::QuantumDot)

Outputs the probability that a given `FockDetection` would have if the `InputType` was `OneParameterInterpolation` with distinguishability `x` for this event. This averages over all possible inputs compatible with the number of lost photons
"""
function p_x_imperfect_source(params_event::SamplingParameters, x, source::QuantumDot) 

    params_event.x = x
    set_parameters!(params_event)
    compute_probability_imperfect_source(params_event, source)

end

"""
    p_x_imperfect_source_update_this_event(event::Event{TIn, TOut},params_event::SamplingParameters, source::QuantumDot) where {TIn<:InputType, TOut <: Union{FockDetection, ThresholdFockDetection}}

Same as above but eating an `Event` so as to keep the previous workflow working.

Still need to create at once the `params_event` at the beginnig of the validation process.
"""
function p_x_imperfect_source_update_this_event(event::Event{TIn, TOut},params_event::SamplingParameters, source::QuantumDot) where {TIn<:InputType, TOut <: Union{FockDetection, ThresholdFockDetection}}

    o = event.output_measurement

    params_event.o = o

    p_x_imperfect_source(params_event, params_event.x, source)

end

p_x_imperfect_source(params_event, 0, source)

ev = params_event_x.ev
ev.output_measurement = FockDetection(ModeOccupation([1,1,1]))
ev

p_x_imperfect_source_update_this_event(ev, params_event_x, source)