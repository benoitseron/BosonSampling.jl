function runs_without_errors(f::Function)

    @test begin
        try
            f()
            true
        catch
            false
        end
    end

end

@testset "Special Matrices" begin

    @testset "rand_gram_matrix_from_orthonormal_basis" begin

        @testset "Basic properties (n=4, r=2)" begin
            n = 4
            r = 2
            S = rand_gram_matrix_from_orthonormal_basis(n, r)

            # Size
            @test size(S) == (n, n)

            # Diagonal elements should be 1 (normalized vectors)
            for i in 1:n
                @test real(S[i,i]) ≈ 1.0 atol=1e-10
                @test imag(S[i,i]) ≈ 0.0 atol=1e-10
            end

            # Hermitian
            @test S ≈ S' atol=1e-10

            # Positive semidefinite (all eigenvalues >= 0)
            eigs = eigvals(S)
            @test all(real(eig) >= -1e-10 for eig in eigs)

            # Rank should be at most r
            @test rank(S) <= r
        end

        @testset "Larger system (n=6, r=3)" begin
            n = 6
            r = 3
            S = rand_gram_matrix_from_orthonormal_basis(n, r)

            @test size(S) == (n, n)
            @test S ≈ S' atol=1e-10
            @test rank(S) <= r
            @test all(real(S[i,i]) ≈ 1.0 for i in 1:n)
        end

        @testset "Rank verification (n=5, r=2)" begin
            n = 5
            r = 2
            S = rand_gram_matrix_from_orthonormal_basis(n, r)

            # Rank should be exactly r (with high probability for random vectors)
            @test rank(S) == r
        end

        @testset "Full-rank case (r >= n)" begin
            # When r >= n the function falls back to a full-rank Gram matrix
            # (rand_gram_matrix), rather than erroring.
            for (n, r) in [(3, 3), (3, 4)]
                S = rand_gram_matrix_from_orthonormal_basis(n, r)
                @test size(S) == (n, n)
                @test S ≈ S' atol=1e-10
                @test rank(S) == n
            end
        end

        @testset "Gram matrix properties" begin
            n = 5
            r = 3
            S = rand_gram_matrix_from_orthonormal_basis(n, r)

            # All off-diagonal elements should have absolute value <= 1
            # (Cauchy-Schwarz inequality)
            for i in 1:n
                for j in 1:n
                    if i != j
                        @test abs(S[i,j]) <= 1.0 + 1e-10
                    end
                end
            end
        end

        @testset "Consistency with rand_gram_matrix_rank" begin
            # Both functions should produce valid Gram matrices
            n = 4
            r = 2

            S1 = rand_gram_matrix_from_orthonormal_basis(n, r)
            S2 = rand_gram_matrix_rank(n, r)

            # Both should have same properties
            @test size(S1) == size(S2) == (n, n)
            @test S1 ≈ S1' atol=1e-10
            @test S2 ≈ S2' atol=1e-10
            @test rank(S1) <= r
            @test rank(S2) <= r
        end

        @testset "Reproducibility" begin
            # Different calls should produce different matrices
            n = 4
            r = 2

            S1 = rand_gram_matrix_from_orthonormal_basis(n, r)
            S2 = rand_gram_matrix_from_orthonormal_basis(n, r)

            # Should be different (with very high probability)
            @test S1 != S2
        end

    end

end
