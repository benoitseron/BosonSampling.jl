@testset "sampling" begin

    @testset "loop" begin
        ### sampling ###

        function loop_tests()

            begin
                n = 3
                m = n

                i = Input{Bosonic}(first_modes(n,m))

                η = 1/sqrt(2) .* ones(m-1)
                η_loss_bs = 0.9 .* ones(m-1)
                η_loss_lines = 0.9 .* ones(m)
                d = Uniform(0, 2pi)
                ϕ = rand(d, m)

            end

            circuit = LossyLoop(m, η, η_loss_bs, η_loss_lines, ϕ).circuit


            p_dark = 0.01
            p_no_count = 0.1

            o = FockSample()
            ev = Event(i,o, circuit)

            BosonSampling.sample!(ev)

            o = DarkCountFockSample(p_dark)
            ev = Event(i,o, circuit)

            BosonSampling.sample!(ev)

            o = RealisticDetectorsFockSample(p_dark, p_no_count)
            ev = Event(i,o, circuit)

            BosonSampling.sample!(ev)


            ###### sample with a new circuit each time ######


            get_sample_loop(LoopSamplingParameters(n = 10 ,input_type = Distinguishable))

            ### method specialisation according to the type of lossy input ###

            n = 6
            m = n

            get_sample_loop(LoopSamplingParameters(n=n, input_type = Distinguishable, η_loss_bs   = nothing, η_loss_lines = 0.9 .* ones(m)))

            get_sample_loop(LoopSamplingParameters(n=n, input_type = Distinguishable, η_loss_bs = 0.9 .* ones(m-1), η_loss_lines = nothing))

            smpl = get_sample_loop(LoopSamplingParameters(n=n, input_type = Distinguishable, η_loss_bs = nothing, η_loss_lines = nothing))

            @test length(smpl.state) == n


        end

        runs_without_errors(loop_tests)
    end

    @testset "closeness of sampling with ideal distribution" begin
        
        n_events = 1000
        n = 2
        m = 2
        interf = Fourier(m)
        TIn = Bosonic
        input_state = Input{TIn}(first_modes(n,m))

        events = generate_experimental_data(n_events = n_events, n = n, m = m, interf = interf, TIn = TIn)

    end


end

n_events = 10000
n = 2
m = 2
interf = Fourier(m)
TIn = Bosonic
input_state = Input{TIn}(first_modes(n,m))

events = generate_experimental_data(n_events = n_events, n = n, m = m, interf = interf, TIn = TIn)

@test tvd_sampled_versus_exact_distribution(events)[1] < 0.05

for n in [2,4,6]

    n_events = 10000
    m = n
    interf = RandHaar(m)
    TIn = Bosonic
    input_state = Input{TIn}(first_modes(n,m))

    events = generate_experimental_data(n_events = n_events, n = n, m = m, interf = interf, TIn = TIn)

    @test tvd_sampled_versus_exact_distribution(events)[1] < 0.05

end

n = 4
n_events = 10000
m = n
interf = RandHaar(m)
TIn = Bosonic
input_state = Input{TIn}(first_modes(n,m))

events = generate_experimental_data(n_events = n_events, n = n, m = m, interf = interf, TIn = TIn)

tvd_sampled_versus_exact_distribution(events)[1] 

