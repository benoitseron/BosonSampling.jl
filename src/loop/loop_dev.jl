include("packages_loop.jl")



@time begin

    n = 9
    m = n

    n_lost = 2

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

    p_x_imperfect_source(params_event, 0, source)

end

