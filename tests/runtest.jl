@testset "permanent ryser" begin

    theoretical_permanent_fourier_matrix = [2,0,-3,0,-5,0,-105,0,81,0,6765,0,175747,0,30375,0,25219857,0,142901109,0,4548104883,0]
    #from "on the permanent of schur matrix", graham

    for n = 1:10
        U = ones(Int64, n, n)
        @test permanent_ryser(U) == factorial(n)
        U = fourier_matrix(n, normalized=false)
        @test permanent_ryser(U) ≈ theoretical_permanent_fourier_matrix[n] atol=1e-2
    end

end

@testset "computing permanent" begin

    for m = 1:10
        U = fourier_matrix(m)
        res = permanent_def(U)

        @test permanent_ryser(U) ≈ res atol=1e-9
        @test multi_dim_ryser(U, gram_matrix_from_x(m, 1)) ≈ abs(res)^2 atol=1e-9
        @test fast_glynn_perm(U) ≈ res atol=1e-9
        @test gurvits(U, 1e-3) ≈ res atol=1e-2

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

            output_events = generate_events(n,n)
            for i = 1:length(output_events)
                if is_forbidden(output_events[i])
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

            output_events = generate_events(n,n)
            for i = 1:length(output_events)
                if is_forbidden(output_events[i])
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

    for n = 3:10

        interf = Fourier(n)
        reflectivity = 1
        distinguishability = 1
        G = GramMatrix{ToyModel}(n, gram_matrix_from_x(n, distinguishability))

        input_clifford_sampler = Input{Bosonic}(first_modes(n,n))
        input_noisy_sampler = Input{ToyModel}(first_modes(n,n), G)

        out_clifford_sampler = cliffords_sampler(input=input_clifford_sampler, interf=interf)
        out_noisy_sampler = noisy_sampling(input=input_noisy_sampler, distinguishability=distinguishability, reflectivity=reflectivity, interf=interf)

        @test !is_forbidden(out_clifford_sampler)
        @test !is_forbidden(out_noisy_sampler)

    end

end
