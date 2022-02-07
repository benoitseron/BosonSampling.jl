using BosonSampling
using Test

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

	@testset "process_probability : HOM" begin
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
		@test proba_partition(U, occupancy_vector) ≈ [0.5,0,0.5] atol = 1e-5

		occupancy_vector = [0, 1]
		@test proba_partition(U, occupancy_vector) ≈ [0.5,0,0.5] atol = 1e-5

		occupancy_vector = [1, 1]
		@test proba_partition(U, occupancy_vector) ≈ [0,0,1.] atol = 1e-5

	end

	@testset "probability distribution of distinguishable photons in modes : HOM" begin

		U = fourier_matrix(2)
		part = [1]
		@test partition_probability_distribution_distinguishable(part, U) ≈ [0.25,0.5,0.25] atol = 1e-8

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

		@test proba_partition_partial(U = U, S = ones(n, n), occupancy_vector = occupancy_vector, input_state = input_state) == proba_partition(U, occupancy_vector, input_state = input_state)
	    @test proba_partition_partial(U = U,  S = Matrix{Float64}(I, n, n), occupancy_vector = occupancy_vector, input_state = input_state) ≈ partition_probability_distribution_distinguishable([1], U; number_photons = n)

	end

end
