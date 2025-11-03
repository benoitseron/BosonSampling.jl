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

@testset "Gram matrix decomposition" begin

    @testset "gram_to_coefficients basic properties" begin
        n = 4
        r = 2
        S = rand_gram_matrix_from_orthonormal_basis(n, r)

        V = gram_to_coefficients(S)

        # Check dimensions
        @test size(V, 1) == n
        @test size(V, 2) <= r
        @test size(V, 2) > 0

        # Reconstruction should match original
        S_reconstructed = reconstruct_gram_matrix(V)
        @test S ≈ S_reconstructed atol=1e-10
    end

    @testset "gram_to_coefficients rank detection" begin
        # Test with different ranks
        for n in [3, 5]
            for r in 1:(n-1)
                S = rand_gram_matrix_from_orthonormal_basis(n, r)
                V = gram_to_coefficients(S)

                # Detected rank should be at least r and at most n
                @test size(V, 2) >= r
                @test size(V, 2) <= n
                @test rank(V) == size(V, 2)

                # Reconstruction should be accurate
                @test S ≈ reconstruct_gram_matrix(V) atol=1e-9
            end
        end
    end

    @testset "gram_to_coefficients with custom atol" begin
        n = 4
        r = 2
        S = rand_gram_matrix_from_orthonormal_basis(n, r)

        # Test with different tolerance
        V = gram_to_coefficients(S, atol=1e-12)
        @test S ≈ reconstruct_gram_matrix(V) atol=1e-10
    end

    @testset "gram_to_coefficients input validation" begin
        n = 3
        r = 2

        # Valid Gram matrix should work
        S_valid = rand_gram_matrix_from_orthonormal_basis(n, r)
        @test_nowarn gram_to_coefficients(S_valid)

        # Invalid Gram matrix (diagonal not 1) should fail
        S_invalid = copy(S_valid)
        S_invalid[1, 1] = 0.5
        @test_throws ArgumentError gram_to_coefficients(S_invalid)

        # Non-Hermitian matrix should fail
        S_non_hermitian = copy(S_valid)
        S_non_hermitian[1, 2] = S_non_hermitian[2, 1] + 0.1
        @test_throws ArgumentError gram_to_coefficients(S_non_hermitian)
    end

    @testset "reconstruct_gram_matrix properties" begin
        n = 5
        r = 3
        S = rand_gram_matrix_from_orthonormal_basis(n, r)
        V = gram_to_coefficients(S)
        S_reconstructed = reconstruct_gram_matrix(V)

        # Reconstructed should be Hermitian
        @test S_reconstructed ≈ S_reconstructed' atol=1e-10

        # Diagonal should be 1
        for i in 1:n
            @test abs(real(S_reconstructed[i, i]) - 1) < 1e-10
            @test abs(imag(S_reconstructed[i, i])) < 1e-10
        end

        # Should be positive semi-definite
        @test minimum(real.(eigvals(S_reconstructed))) >= -1e-10
    end

end 

