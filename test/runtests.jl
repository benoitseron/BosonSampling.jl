using BosonSampling
using Plots
using ProgressMeter
using Distributions
using Random
using Test
using ArgCheck
using StatsBase
using ColorSchemes
using Interpolations
using Dierckx
using LinearAlgebra
using PrettyTables
using LaTeXStrings
using JLD
using AutoHashEquals
using LinearRegression
using DataStructures

include("tools.jl")

@testset "BosonSampling.jl" begin

    @testset "scattering matrix" begin

        input_state = [2,0]
	    output_state = [1,1]
	    U = matrix_test(2)
	    @test scattering_matrix(U, input_state, output_state) == [1 2; 1 2]

	    input_state = [1,1]
	    output_state = [0,2]
	    U = matrix_test(2)
	    @test scattering_matrix(U, input_state, output_state) == [2 2; 4 4]

	end

	@testset "process_probability: HOM" begin

		@test process_probability(fourier_matrix(2), [1,1], [0,2]) ≈ 0.5 atol = 1e-8
	    @test process_probability(fourier_matrix(2), [1,1], [2,0]) ≈ 0.5 atol = 1e-8
	    @test process_probability(fourier_matrix(2), [1,1], [1,1]) ≈ 0. atol = 1e-8

		S_bosonic = ones(2,2)
		S_dist = [1 0; 0 1]

		@test process_probability_partial(fourier_matrix(2), S_bosonic, [1,1], [0,2]) ≈ 0.5 atol = 1e-8
		@test process_probability_partial(fourier_matrix(2), S_bosonic, [1,1], [2,0]) ≈ 0.5 atol = 1e-8
		@test process_probability_partial(fourier_matrix(2), S_bosonic, [1,1], [1,1]) ≈ 0. atol = 1e-8
		@test process_probability_partial(fourier_matrix(2), S_dist, [1,1], [0,2]) ≈ 0.25 atol = 1e-8
		@test process_probability_partial(fourier_matrix(2), S_dist, [1,1], [2,0]) ≈ 0.25 atol = 1e-8
		@test process_probability_partial(fourier_matrix(2), S_dist, [1,1], [1,1]) ≈ 0.5 atol = 1e-8

	end

	@testset "probability distribution of photons in modes : HOM" begin

		U = fourier_matrix(2)

		occupancy_vector = [1, 0]
		@test proba_partition_bosonic(U = U, occupancy_vector = occupancy_vector) ≈ [0.5,0,0.5] atol = 1e-5

		occupancy_vector = [0, 1]
		@test proba_partition_bosonic(U = U, occupancy_vector = occupancy_vector) ≈ [0.5,0,0.5] atol = 1e-5

		occupancy_vector = [1, 1]
		@test proba_partition_bosonic(U = U, occupancy_vector = occupancy_vector) ≈ [0,0,1.] atol = 1e-5

	end

	@testset "probability distribution of distinguishable photons in modes : HOM" begin

	    U = fourier_matrix(2)
	    part = [1,0]
	    @test proba_partition_distinguishable(occupancy_vector = part, U = U) ≈ [0.25,0.5,0.25] atol = 1e-8

	end

	@testset "partial distinguishability partitions" begin

		m = 10
		n = 4
		input_state = zeros(Int, m)
		occupancy_vector = zeros(Int, m)
		for i=1:n
			input_state[i] = 1
		end
		occupancy_vector[1] = 1

		U = fourier_matrix(10)

		@test proba_partition_partial(U = U, S = ones(n, n), occupancy_vector = occupancy_vector, input_state = input_state) == proba_partition_bosonic(U = U, occupancy_vector = occupancy_vector, input_state = input_state)
	    @test proba_partition_partial(U = U,  S = Matrix{Float64}(I, n, n), occupancy_vector = occupancy_vector, input_state = input_state) ≈ proba_partition_distinguishable(occupancy_vector = occupancy_vector, U = U, input_state = input_state)

	end

	@testset "theoretical_distribution" begin

		n = 3
		occupation = random_mode_occupation_collisionless(n,n)
		i = Input{Bosonic}(occupation)
		interf = Fourier(i.r.m)
		p_theo = theoretical_distribution(input=i, interf=interf)

		@test sum(p_theo) ≈ 1 atol=1e-9
		@test all(p -> p>=0, p_theo)

		output_events = output_mode_occupation(n,n)
	    for i = 1:length(output_events)
	        if check_suppression_law(output_events[i])
	            @test p_theo[i] ≈ 0 atol=1e-6
	        end
	    end

	end

	@testset "noisy statistics" begin

		n = 3
		occupation = random_mode_occupation_collisionless(n,n)
    	i = Input{OneParameterInterpolation}(occupation, 0.5)
	    interf = Fourier(i.r.m)
	    res = noisy_distribution(input=i, loss=0.5, interf=interf)

		p_exact = res[1]
		p_approx = res[2]
		p_samp = res[3]

		@test sum(p_exact) ≈ 1 atol=1e-9
		@test sum(p_approx) ≈ 1 atol=1e-9
		@test sum(p_samp) ≈ 1 atol=1e-9
		@test all(p->p>=0, p_exact)
		@test all(p->p>=0, p_approx)
		@test all(p->p>=0, p_samp)

		i = Input{Bosonic}(occupation)
		res = noisy_distribution(input=i, loss=0.999, interf=interf, approx=false, samp=false)
		p_exact = res[1]
		output_events = output_mode_occupation(n,n)
		for i = 1:length(output_events)
            if check_suppression_law(output_events[i])
		        @test p_exact[i] ≈ 0 atol=1e-6
		    end
		end

	end

	@testset "suppression law boson samplers" begin

		n = 3
	    interf = Fourier(n)

		input_clifford_sampler = Input{Bosonic}(first_modes(n,n))
	   	input_noisy_sampler = Input{OneParameterInterpolation}(first_modes(n,n), 1.0)

	    out_clifford_sampler = cliffords_sampler(input=input_clifford_sampler, interf=interf)
		out_noisy_sampler = noisy_sampler(input=input_noisy_sampler, loss=1.0, interf=interf)

	    @test !check_suppression_law(out_clifford_sampler)
		@test !check_suppression_law(out_noisy_sampler)

	end

	@testset "check analytical counter example interferometer" begin

		n = 7
		r = 2

		complete_interferometer = Matrix{ComplexF16}(I,n,n)
		bottom_dft = Matrix{ComplexF64}(I,n,n)
		bottom_dft[r+1:n, r+1:n] = fourier_matrix(n-r)
		complete_interferometer *= bottom_dft

		for top_mode in 1:r
			complete_interferometer *= beam_splitter_modes(in_up=top_mode, in_down=top_mode+r, out_up=top_mode, out_down=top_mode+r,transmission_amplitude=sqrt(r/n), n=n)
		end

		complete_interferometer = transpose(complete_interferometer)

		circuit = Circuit(n)
		add_element!(circuit=circuit, interf=Fourier(n-r), target_modes=[i for i in r+1:n])
		add_element!(circuit=circuit, interf=BeamSplitter(sqrt(r/n)), target_modes=[3,1])
		add_element!(circuit=circuit, interf=BeamSplitter(sqrt(r/n)), target_modes=[4,2])

		@test complete_interferometer == circuit.U

	end


	@testset "examples usage" begin
		@test include("example_usage.jl")
	end

	include("circuits_and_loss.jl")
	include("partitions.jl")
	include("sources.jl")
	include("threshold.jl")

end
