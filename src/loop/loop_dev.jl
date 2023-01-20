include("packages_loop.jl")

using TimerOutputs
using Profile
using ProfileView


func() = begin

    n = 8
    m = n

    n_lost = 1

    source = QuantumDot(13.5 / 80)

    d = Uniform(0,2pi)
    ϕ = nothing # rand(d,m)

    η_loss_lines = 0.9 * ones(m)
    η_loss_bs =  0.9 * ones(m-1)

    params = LoopSamplingParameters(n=n, η = η_thermalization(n), η_loss_bs = η_loss_bs, η_loss_lines = η_loss_lines, ϕ = ϕ)

    params_event = convert(SamplingParameters, params)
    params_event.o = ThresholdFockDetection(ThresholdModeOccupation(first_modes(n-n_lost, m).state))

    set_parameters!(params_event)

    ev = params_event.ev

    # can also just use 

    compute_probability!(params_event) 

    ### with the validation formalism now ###

    @timeit to "p_x_imperfect_source" p_x_imperfect_source(params_event, 0, source)

end


n = 10
m = n
n_lost = 1

d = Uniform(0,2pi)
ϕ = nothing # rand(d,m)
η_loss_lines = nothing # 0.2 * ones(m)
η_loss_bs = nothing #0.2 * ones(m-1)

params = LoopSamplingParameters(n=n, η = η_thermalization(n), η_loss_bs = η_loss_bs, η_loss_lines = η_loss_lines, ϕ = ϕ)

params_event = convert(SamplingParameters, params)
params_event.o = ThresholdFockDetection(ThresholdModeOccupation(first_modes(n-n_lost, m).state))

set_parameters!(params_event)

ev = params_event.ev

@time compute_threshold_detection_probability(ev)
@time compute_probability!(ev)

params_event.o = FockDetection(ModeOccupation(ones(n)))

set_parameters!(params_event)

ev = params_event.ev

@time compute_probability!(ev)

@time ryser(ev.interferometer.U[1:10,1:10])

