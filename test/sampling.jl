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

        end

        runs_without_errors(loop_tests)
    end

end
