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
end

@testset "theoretical distribution" begin

    @testset "normalization" begin

        for n = 2:5
            occupation = ModeOccupation(random_occupancy(n, n^2))
            input = Input{Bosonic}(occupation)
            interf = RandHaar(input.r.m)
            p_theo = theoretical_distribution(input=input, distinguishability=1, interf=interf, gram_matrix=input.G)

            @test sum(p_theo) ≈ 1 atol=1e-9
        end

    end

    @testset "positivity" begin

        p = Progress(4, 1, "check positivity...", 50)
        for n = 2:5
            occupation = ModeOccupation(random_occupancy(n, n^2))
            input = Input{Bosonic}(occupation)
            interf = RandHaar(input.r.m)
            p_theo = theoretical_distribution(input=input, distinguishability=1, interf=interf, gram_matrix=input.G)

            @test all(p->p>=0, p_theo)
            next!(p)
        end

    end

    @testset "suppression law" begin

        for n = 2:5
            occupation = ModeOccupation(random_occupancy(n,n))
            input = Input{Bosonic}(occupation)
            interf = Fourier(input.r.n)
            p_theo = theoretical_distribution(input=input, distinguishability=1, interf=interf, gram_matrix=input.G)

            output_events = output_mode_occupation(n,n)
            for i = 1:length(output_events)
                if check_suppression_law(output_events[i])
                    @test p_theo[i] ≈ 0 atol=1e-6
                end
            end
        end

    end
end

@testset "noisy statistics" begin

    @testset "normalization" begin

        for n = 2:4
            occupation = ModeOccupation(random_occupancy(n, 2n))
            input = Input{RandomModel}(occupation)
            interf = RandHaar(input.r.m)
            p_exact, p_approx, p_samp = noisy_distribution(input=input, distinguishability=0.5, reflectivity=0.5, interf=interf)

            @test sum(p_exact) ≈ 1 atol=1e-9
            @test sum(p_approx) ≈ 1 atol=1e-9
            @test sum(p_samp) ≈ 1 atol=1e-9
        end

    end

    @testset "positivity" begin

        for n = 2:4
            occupation = ModeOccupation(random_occupancy(n, 2n))
            input = Input{RandomModel}(occupation)
            interf = RandHaar(input.r.m)
            p_exact, p_approx, p_samp = noisy_distribution(input=input, distinguishability=0.5, reflectivity=0.5, interf=interf)

            @test all(p->p>=0, p_exact)
            @test all(p->p>=0, p_approx)
            @test all(p->p>=0, p_samp)
        end

    end

    @testset "suppression law" begin

        # https://arxiv.org/pdf/1002.5038.pdf

        for n = 2:4

            occupation = ModeOccupation(random_occupancy(n,n))
            input = Input{Bosonic}(occupation)
            interf = Fourier(input.r.m)
            O = noisy_distribution(input=input, distinguishability=1, reflectivity=0.5, interf=interf, approx=false, samp=false)
            p_exact = O[1]

            output_events = output_mode_occupation(n,n)
            for i = 1:length(output_events)
                if check_suppression_law(output_events[i])
                    @test p_exact[i] ≈ 0 atol=1e-6
                end
            end

        end

    end

end

@testset "consistency" begin

    for n = 2:4
        occupation = ModeOccupation(random_occupancy(n, 2n))
        input = Input{Bosonic}(occupation)

        for i = 1:10
            interf = RandHaar(input.r.m)
            p_theo = theoretical_distribution(input=input, distinguishability=1, interf=interf, gram_matrix=input.G)
            O = noisy_distribution(input=input, distinguishability=1, reflectivity=0.999, interf=interf, approx=false, samp=false)
            p_exact = O[1]

            for j = 1:length(p_theo)
                @test abs(p_theo[i]-p_exact[i]) ≈ 0 atol=1e-9
            end
        end
    end

end

@testset "suppression law boson samplers" begin
	@warn "need to include Permanents.jl to export fast_glynn_pern"

    for n = 3:10

        interf = Fourier(n)
        reflectivity = 1
        distinguishability = 1
        G = GramMatrix{ToyModel}(n, gram_matrix_toy_model(n, distinguishability))

        input_clifford_sampler = Input{Bosonic}(first_modes(n,n))
        input_noisy_sampler = Input{ToyModel}(first_modes(n,n), G)

        out_clifford_sampler = cliffords_sampler(input=input_clifford_sampler, interf=interf)
        out_noisy_sampler = noisy_sampling(input=input_noisy_sampler, distinguishability=distinguishability, reflectivity=reflectivity, interf=interf)

        @test !check_suppression_law(out_clifford_sampler)
        @test !check_suppression_law(out_noisy_sampler)

    end

end