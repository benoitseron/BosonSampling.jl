"""
Test sample_multiple function
"""

@testset "sample_multiple function" begin

    # Setup
    n = 3
    m = 5
    n_samples = 50

    @testset "UserDefinedGramMatrix" begin
        S = rand_gram_matrix_from_orthonormal_basis(n, 2)
        input = Input{UserDefinedGramMatrix}(first_modes(n, m), S)
        interf = RandHaar(m)

        samples = sample_multiple(input, interf, n_samples)

        @test length(samples) == n_samples
        @test all(length(s) == m for s in samples)
        @test all(sum(s) == n for s in samples)
        @test all(all(x >= 0 for x in s) for s in samples)  # Non-negative counts
    end

    @testset "OneParameterInterpolation" begin
        input = Input{OneParameterInterpolation}(first_modes(n, m), 0.7)
        interf = RandHaar(m)

        samples = sample_multiple(input, interf, n_samples)

        @test length(samples) == n_samples
        @test all(length(s) == m for s in samples)
        @test all(sum(s) == n for s in samples)
    end

    @testset "Bosonic (fallback)" begin
        input = Input{Bosonic}(first_modes(n, m))
        interf = RandHaar(m)

        samples = sample_multiple(input, interf, n_samples)

        @test length(samples) == n_samples
        @test all(length(s) == m for s in samples)
        @test all(sum(s) == n for s in samples)
    end

    @testset "Distinguishable (fallback)" begin
        input = Input{Distinguishable}(first_modes(n, m))
        interf = RandHaar(m)

        samples = sample_multiple(input, interf, n_samples)

        @test length(samples) == n_samples
        @test all(length(s) == m for s in samples)
        @test all(sum(s) == n for s in samples)
    end

    @testset "Edge cases" begin
        # Single sample
        input = Input{UserDefinedGramMatrix}(first_modes(n, m), rand_gram_matrix_from_orthonormal_basis(n, 2))
        interf = RandHaar(m)

        samples = sample_multiple(input, interf, 1)
        @test length(samples) == 1
        @test sum(samples[1]) == n

        # Many samples (verify consistency)
        samples_many = sample_multiple(input, interf, 200)
        @test length(samples_many) == 200
        @test all(sum(s) == n for s in samples_many)
    end

    @testset "Consistency with sample!" begin
        # Verify sample_multiple gives same distribution as repeated sample!
        input = Input{OneParameterInterpolation}(first_modes(n, m), 0.5)
        interf = Fourier(m)  # Use deterministic interferometer

        # This test just verifies both methods work, not that they're identical
        # (since sampling is random)
        samples_multiple = sample_multiple(input, interf, 10)

        samples_single = []
        for i in 1:10
            ev = Event(input, FockSample(), interf)
            sample!(ev)
            push!(samples_single, ev.output_measurement.s.state)
        end

        @test all(sum(s) == n for s in samples_multiple)
        @test all(sum(s) == n for s in samples_single)
    end

end
