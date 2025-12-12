"""
Tests for the Clifford sampler to verify it produces distributions
consistent with compute_probability! for all cases including bunching.

These tests verify the fix for the matrix convention issue where the
standard Clifford algorithm uses U[output, input] but this package
uses U[input, output] convention.
"""

using Test
using BosonSampling
using Random
using StatsBase
using LinearAlgebra

@testset "Clifford Sampler" begin

    @testset "HOM dip (collision-free baseline)" begin
        # Standard HOM test - collision-free, so both conventions work
        U = fourier_matrix(2)
        @test process_probability(U, [1,1], [0,2]) ≈ 0.5 atol = 1e-8
        @test process_probability(U, [1,1], [2,0]) ≈ 0.5 atol = 1e-8
        @test process_probability(U, [1,1], [1,1]) ≈ 0. atol = 1e-8
    end

    @testset "Sampler vs compute_probability (2 photons, 4 modes)" begin
        Random.seed!(42)
        n, m = 2, 4
        interf = RandHaar(m)
        input = Input{Bosonic}(first_modes(n, m))

        # Get unique outputs
        mode_lists = output_mode_occupation(n, m)
        output_events = unique([mode_occupancy_to_occupancy_vector(ml, m) for ml in mode_lists])

        # Compute exact probabilities via compute_probability!
        exact_probs = Dict{Vector{Int}, Float64}()
        for out_state in output_events
            mode_occ = ModeOccupation(out_state)
            output = FockDetection(mode_occ)
            ev = Event(input, output, interf)
            compute_probability!(ev)
            exact_probs[out_state] = real(ev.proba_params.probability)
        end

        # Sample many times
        n_samples = 20000
        sample_counts = Dict{Vector{Int}, Int}()
        for _ in 1:n_samples
            s = cliffords_sampler(input=input, interf=interf)
            sample_counts[s] = get(sample_counts, s, 0) + 1
        end

        # Compute TVD
        tvd = 0.0
        for out_state in output_events
            p_exact = exact_probs[out_state]
            p_samp = get(sample_counts, out_state, 0) / n_samples
            tvd += abs(p_exact - p_samp)
        end
        tvd /= 2

        @test tvd < 0.02
    end

    @testset "Sampler vs compute_probability (3 photons, 5 modes)" begin
        Random.seed!(123)
        n, m = 3, 5
        interf = RandHaar(m)
        input = Input{Bosonic}(first_modes(n, m))

        # Get unique outputs
        mode_lists = output_mode_occupation(n, m)
        output_events = unique([mode_occupancy_to_occupancy_vector(ml, m) for ml in mode_lists])

        # Compute exact probabilities
        exact_probs = Dict{Vector{Int}, Float64}()
        for out_state in output_events
            mode_occ = ModeOccupation(out_state)
            output = FockDetection(mode_occ)
            ev = Event(input, output, interf)
            compute_probability!(ev)
            exact_probs[out_state] = real(ev.proba_params.probability)
        end

        # Sample many times
        n_samples = 30000
        sample_counts = Dict{Vector{Int}, Int}()
        for _ in 1:n_samples
            s = cliffords_sampler(input=input, interf=interf)
            sample_counts[s] = get(sample_counts, s, 0) + 1
        end

        # Compute TVD
        tvd = 0.0
        for out_state in output_events
            p_exact = exact_probs[out_state]
            p_samp = get(sample_counts, out_state, 0) / n_samples
            tvd += abs(p_exact - p_samp)
        end
        tvd /= 2

        @test tvd < 0.03
    end

    @testset "Bunching outputs (collision case)" begin
        # This test specifically checks bunching outputs where the
        # matrix convention matters most
        Random.seed!(999)
        n, m = 2, 3
        interf = RandHaar(m)
        input = Input{Bosonic}(first_modes(n, m))

        # Test all bunching outputs [2,0,0], [0,2,0], [0,0,2]
        bunching_outputs = [[2,0,0], [0,2,0], [0,0,2]]

        for bunched_output in bunching_outputs
            # Compute exact probability
            output = FockDetection(ModeOccupation(bunched_output))
            ev = Event(input, output, interf)
            compute_probability!(ev)
            p_exact = real(ev.proba_params.probability)

            # Sample and count
            n_samples = 50000
            count_bunched = 0
            for _ in 1:n_samples
                s = cliffords_sampler(input=input, interf=interf)
                if s == bunched_output
                    count_bunched += 1
                end
            end
            p_sampled = count_bunched / n_samples

            @test abs(p_exact - p_sampled) < 0.01
        end
    end

    @testset "Various configurations TVD" begin
        # Test multiple (n, m) configurations to ensure robustness
        configs = [(2, 2), (2, 4), (2, 6), (3, 3), (3, 4), (4, 4)]

        for (n, m) in configs
            Random.seed!(42 + n*10 + m)
            interf = RandHaar(m)
            input = Input{Bosonic}(first_modes(n, m))

            # Get unique outputs
            mode_lists = output_mode_occupation(n, m)
            output_events = unique([mode_occupancy_to_occupancy_vector(ml, m) for ml in mode_lists])

            # Compute exact probabilities
            exact_probs = Dict{Vector{Int}, Float64}()
            for out_state in output_events
                mode_occ = ModeOccupation(out_state)
                output = FockDetection(mode_occ)
                ev = Event(input, output, interf)
                compute_probability!(ev)
                exact_probs[out_state] = real(ev.proba_params.probability)
            end

            # Sample
            n_samples = 15000
            sample_counts = Dict{Vector{Int}, Int}()
            for _ in 1:n_samples
                s = cliffords_sampler(input=input, interf=interf)
                sample_counts[s] = get(sample_counts, s, 0) + 1
            end

            # Compute TVD
            tvd = 0.0
            for out_state in output_events
                p_exact = exact_probs[out_state]
                p_samp = get(sample_counts, out_state, 0) / n_samples
                tvd += abs(p_exact - p_samp)
            end
            tvd /= 2

            @test tvd < 0.03
        end
    end

    @testset "Sample validity" begin
        # Verify all samples have correct photon number and valid modes
        Random.seed!(12345)
        n, m = 3, 6
        input = Input{Bosonic}(first_modes(n, m))
        interf = RandHaar(m)

        for _ in 1:100
            s = cliffords_sampler(input=input, interf=interf)

            # Check photon conservation
            @test sum(s) == n

            # Check valid mode occupation
            @test length(s) == m
            @test all(s .>= 0)
        end
    end

    @testset "Deterministic with seed" begin
        # Verify reproducibility with same seed
        Random.seed!(42)
        n, m = 2, 4
        interf = RandHaar(m)
        input = Input{Bosonic}(first_modes(n, m))

        Random.seed!(100)
        s1 = cliffords_sampler(input=input, interf=interf)

        Random.seed!(100)
        s2 = cliffords_sampler(input=input, interf=interf)

        @test s1 == s2
    end

end
