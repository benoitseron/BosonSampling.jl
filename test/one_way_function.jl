# ═══════════════════════════════════════════════════════════════════════
# One-way function via boson sampling — tests
# ═══════════════════════════════════════════════════════════════════════
# The estimator lives in src/boson_samplers/one_way_function.jl and is compiled
# into BosonSampling, so no include of the source is needed here.  BosonSampling auto-exports
# every module identifier (see the loop at the end of src/BosonSampling.jl), so
# `using BosonSampling` brings the whole API — including internals such as
# `_harmonic`, `_compute_N`, `_glynn_fused!` — into scope.
using LinearAlgebra
using Random
using Statistics: mean, std
using BosonSampling
using Permanents: ryser
using Test

# ═══════════════════════════════════════════════════════════════════════
# Helper: enumerate ALL output configs (including collisions)
# ═══════════════════════════════════════════════════════════════════════
function enumerate_all_outputs(m::Int, n::Int)
    results = Vector{Vector{Int}}()
    _enum_all!(results, Int[], m, n)
    return results
end

function _enum_all!(results, current, m, n_remaining)
    if length(current) == m
        if n_remaining == 0
            push!(results, copy(current))
        end
        return
    end
    for s in 0:n_remaining
        push!(current, s)
        _enum_all!(results, current, m, n_remaining - s)
        pop!(current)
    end
end

# Exact boson sampling probability Pr[output] for input |1,1,...,1,0,...,0⟩
# Formula: Pr[s] = |Per(U_sub)|^2 / (s_1! ... s_m!)
# where U_sub is n×n: rows = output modes (repeated), cols = input modes 1:n
function exact_probability(U::AbstractMatrix, n::Int, output::Vector{Int})
    m = length(output)
    out_cols = Int[]
    for j in 1:m
        for _ in 1:output[j]
            push!(out_cols, j)
        end
    end
    U_sub = U[out_cols, 1:n]
    perm_val = ryser(U_sub)
    denom = prod(factorial(s) for s in output)
    return abs2(perm_val) / denom
end

# Exact characteristic function by definition (Eq. 6):
# φ(k) = sum_s Pr[s] exp(-2πik f(s)/N)
function exact_char_func(U::AbstractMatrix, n::Int, m::Int, omega, N, k)
    all_out = enumerate_all_outputs(m, n)
    result = complex(0.0)
    for s in all_out
        p = exact_probability(U, n, s)
        fs = sum(omega[i] * s[i] for i in 1:m)
        result += p * exp(-2π * im * k * fs / N)
    end
    return result
end

# Exact CDF S(x0) = sum_{s: f(s) <= x0} Pr[s]
function exact_cdf(U::AbstractMatrix, n::Int, m::Int, omega, x0)
    all_out = enumerate_all_outputs(m, n)
    result = 0.0
    for s in all_out
        fs = sum(omega[i] * s[i] for i in 1:m)
        if fs <= x0
            result += exact_probability(U, n, s)
        end
    end
    return result
end


@testset "one-way function (boson sampling OWF)" begin

# ═══════════════════════════════════════════════════════════════════════
# TEST 1: Glynn single-sample estimator
# ═══════════════════════════════════════════════════════════════════════
@testset "Glynn estimator" begin

    # 1a. Exhaustive average over all x ∈ {-1,1}^n must equal Per(A)
    #     E_{x}[Gly_x(A)] = Per(A)  (Eq. 21 of Aaronson-Hance)
    @testset "Exhaustive average = Per(A)" begin
        for trial in 1:5
            Random.seed!(trial)
            n = rand(2:5)
            A = randn(ComplexF64, n, n)

            # Exhaustive: enumerate all x ∈ {-1,1}^n
            avg = complex(0.0)
            for bits in 0:(2^n - 1)
                x = [2 * ((bits >> i) & 1) - 1 for i in 0:n-1]
                avg += glynn_single(A, x)
            end
            avg /= 2^n

            perm_exact = ryser(A)
            @test abs(avg - perm_exact) < 1e-10 * max(1.0, abs(perm_exact))
        end
    end

    # 1b. Known permanents
    @testset "Known permanents" begin
        # Identity: Per(I_n) = 1
        for n in 2:5
            I_n = Matrix{ComplexF64}(I, n, n)
            avg = complex(0.0)
            for bits in 0:(2^n - 1)
                x = [2 * ((bits >> i) & 1) - 1 for i in 0:n-1]
                avg += glynn_single(I_n, x)
            end
            avg /= 2^n
            @test abs(avg - 1.0) < 1e-12
        end

        # Ones matrix: Per(J_n) = n!
        for n in 2:5
            J = ones(ComplexF64, n, n)
            avg = complex(0.0)
            for bits in 0:(2^n - 1)
                x = [2 * ((bits >> i) & 1) - 1 for i in 0:n-1]
                avg += glynn_single(J, x)
            end
            avg /= 2^n
            @test abs(avg - factorial(n)) < 1e-10
        end
    end

    # 1c. Bound: |Gly_x(A)| ≤ ||A||^n (Proposition 1 of Aaronson-Hance)
    @testset "Operator norm bound" begin
        Random.seed!(42)
        for _ in 1:20
            n = rand(2:6)
            A = randn(ComplexF64, n, n)
            x = rand([-1, 1], n)
            norm_A = opnorm(A)
            @test abs(glynn_single(A, x)) <= norm_A^n + 1e-10
        end
    end

    # 1d. For unitary matrices |Gly_x(U)| ≤ 1 (since ||U|| = 1)
    @testset "Unitary bound |Gly_x(U)| ≤ 1" begin
        Random.seed!(42)
        for _ in 1:20
            n = rand(2:6)
            U = RandHaar(n).U
            x = rand([-1, 1], n)
            @test abs(glynn_single(U, x)) <= 1.0 + 1e-10
        end
    end

    # 1e. Statistical convergence: Monte Carlo average → Per(A) as M → ∞
    @testset "Monte Carlo convergence" begin
        Random.seed!(42)
        n = 4
        A = randn(ComplexF64, n, n)
        perm_exact = ryser(A)

        M = 100_000
        mc_avg = complex(0.0)
        for _ in 1:M
            x = rand([-1, 1], n)
            mc_avg += glynn_single(A, x)
        end
        mc_avg /= M

        # With M=100k, error should be ~||A||^n / sqrt(M)
        rel_tol = 5 * opnorm(A)^n / sqrt(M)
        @test abs(mc_avg - perm_exact) < rel_tol
    end
end


# ═══════════════════════════════════════════════════════════════════════
# TEST 2: Importance sampling distribution q(k) (Proposition 2)
# ═══════════════════════════════════════════════════════════════════════
@testset "Importance sampling distribution" begin

    # 2a. t(k) matches Eq. 17 for specific values
    @testset "t(k) definition (Eq. 17)" begin
        N = 20
        @test t_weight(0, N) == 1.0
        @test t_weight(1, N) == 1.0 / 1
        @test t_weight(5, N) == 1.0 / 5
        @test t_weight(9, N) == 1.0 / 9    # N/2 - 1 = 9
        @test t_weight(10, N) == 1.0 / 10   # N/2 = 10 → 1/(N-k) = 1/10
        @test t_weight(15, N) == 1.0 / 5    # 1/(N-15) = 1/5
        @test t_weight(19, N) == 1.0 / 1    # 1/(N-19) = 1/1
    end

    # 2b. q(k) sums to 1
    @testset "q(k) is a valid PMF" begin
        for N in [10, 17, 50, 100]
            𝒩 = normalization_constant(N)
            total = sum(t_weight(k, N) / 𝒩 for k in 0:N-1)
            @test abs(total - 1.0) < 1e-12
        end
    end

    # 2c. _harmonic matches exact sum for moderate k
    @testset "_harmonic accuracy" begin
        for k in [1, 10, 100, 1000, 5000]
            exact = sum(1.0 / i for i in 1:k)
            @test abs(_harmonic(k) - exact) < 1e-12
        end
        # Asymptotic branch: verify relative accuracy for large k
        for k in [20_000, 100_000, 1_000_000]
            h = _harmonic(k)
            # H_k ~ ln(k) + γ, so it should be positive and growing
            @test h > log(k)
            @test h < log(k) + 1.0
        end
    end

    # 2d. normalization_constant matches direct sum for small N
    @testset "normalization_constant matches direct sum" begin
        for N in [3, 10, 17, 50, 100, 500]
            direct = sum(t_weight(k, N) for k in 0:N-1)
            @test abs(normalization_constant(N) - direct) < 1e-10
        end
    end

    # 2e. t(k) ~ 1/min(k, N-k)
    @testset "t(k) ∝ 1/min(k, N-k)" begin
        for N in [16, 17, 50, 51]
            for k in 1:N-1
                mk = min(k, N - k)
                @test t_weight(k, N) <= 1.0 / mk + 1e-14
                @test t_weight(k, N) >= 1.0 / (mk + 1) - 1e-14
            end
        end
    end

    # 2f. _sample_reciprocal: empirical distribution matches ∝ 1/k
    @testset "_sample_reciprocal distribution" begin
        Random.seed!(42)
        K = 20
        M = 200_000
        counts = zeros(Int, K)
        for _ in 1:M
            k = _sample_reciprocal(K)
            counts[k] += 1
        end

        H_K = sum(1.0 / k for k in 1:K)
        for k in 1:K
            expected_freq = M / (k * H_K)
            if expected_freq > 100
                @test abs(counts[k] - expected_freq) / expected_freq < 0.1
            end
        end
    end

    # 2g. sample_fourier_mode: empirical distribution matches q(k)
    @testset "sample_fourier_mode distribution" begin
        Random.seed!(42)
        N = 20
        M = 300_000
        counts = zeros(Int, N)
        for _ in 1:M
            k = sample_fourier_mode(N)
            counts[k + 1] += 1  # 0-indexed → 1-indexed
        end

        𝒩 = normalization_constant(N)
        for k in 0:N-1
            expected_freq = M * t_weight(k, N) / 𝒩
            if expected_freq > 100
                rel_err = abs(counts[k + 1] - expected_freq) / expected_freq
                @test rel_err < 0.1
            end
        end
    end
end


# ═══════════════════════════════════════════════════════════════════════
# TEST 3: Geometric sum G_N(k; x0) (Eq. 9)
# ═══════════════════════════════════════════════════════════════════════
@testset "Geometric sum G_N" begin

    # 3a. Direct summation check: G_N(k, x0) = sum_{j=0}^{x0} exp(2πijk/N)
    @testset "Matches direct summation" begin
        for N in [10, 17, 50]
            for k in 0:N-1
                for x0 in [0, 1, N ÷ 3, N ÷ 2, N - 1]
                    direct = sum(exp(2π * im * j * k / N) for j in 0:x0)
                    @test abs(G_N(k, x0, N) - direct) < 1e-10
                end
            end
        end
    end

    # 3b. k=0: G_N(0, x0) = x0 + 1
    @testset "k=0 case" begin
        for N in [10, 50]
            for x0 in [0, 5, N - 1]
                @test abs(G_N(0, x0, N) - (x0 + 1)) < 1e-12
            end
        end
    end

    # 3c. x0 = N-1 and k ≠ 0: sum of all N-th roots of unity = 0
    @testset "Full sum = 0 for k ≠ 0" begin
        for N in [10, 17, 50]
            for k in 1:N-1
                @test abs(G_N(k, N - 1, N)) < 1e-10
            end
        end
    end

    # 3d. Bound |G_N(k, x0)| ≤ N / (2 min(k, N-k)) for k ≠ 0 (Eq. 25)
    @testset "Bound (Eq. 25)" begin
        N = 50
        for k in 1:N-1
            for x0 in [0, N ÷ 4, N ÷ 2, N - 1]
                bound = N / (2 * min(k, N - k))
                @test abs(G_N(k, x0, N)) <= bound + 1e-10
            end
        end
    end
end


# ═══════════════════════════════════════════════════════════════════════
# TEST 4: Diagonal phase matrix D_k
# ═══════════════════════════════════════════════════════════════════════
@testset "Diagonal phase matrix D_k" begin

    omega = [1.0, 3.0, 9.0, 27.0]  # base-3 (n+1=3 for n=2)
    N = 55

    # 4a. D_0 = Identity
    @testset "D_0 = I" begin
        D0 = diag_Dk(0, omega, N)
        @test D0 ≈ I(4)
    end

    # 4b. Diagonal entries match definition
    @testset "Diagonal entries" begin
        k = 3
        Dk = diag_Dk(k, omega, N)
        for j in 1:4
            expected = exp(-2π * im * k * omega[j] / N)
            @test abs(Dk[j, j] - expected) < 1e-14
        end
    end

    # 4c. D_k is unitary (diagonal with |entries| = 1)
    @testset "D_k is unitary" begin
        for k in 0:N-1
            Dk = diag_Dk(k, omega, N)
            @test norm(Dk' * Dk - I(4)) < 1e-12
        end
    end

    # 4d. U† D_k U is unitary when U is unitary
    @testset "U† D_k U is unitary" begin
        U = RandHaar(4).U
        for k in 0:N-1
            Dk = diag_Dk(k, omega, N)
            A = U' * Dk * U
            @test norm(A' * A - I(4)) < 1e-12
        end
    end
end


# ═══════════════════════════════════════════════════════════════════════
# TEST 5: Characteristic function identity (THE KEY THEOREM)
# φ(k) = Per(U_in† D_k U_in) = E_{s~D_U}[exp(-2πik f(s)/N)]  (Eq. 12)
# Now using base-(n+1) weights for collision-full regime.
# ═══════════════════════════════════════════════════════════════════════
@testset "Characteristic function identity (Eq. 12)" begin

    for (n, m) in [(1, 3), (1, 5), (2, 3), (2, 4), (2, 5), (3, 4), (3, 5)]
        @testset "n=$n, m=$m" begin
            Random.seed!(100 + n + 10m)
            U = RandHaar(m).U
            base = n + 1
            omega = [Float64(base^(i - 1)) for i in 1:m]
            N = n * base^(m - 1) + 1

            for k in 0:min(N - 1, 20)
                # LHS: Per(U_in† D_k U_in) via exact permanent
                Dk = diag_Dk(k, omega, N)
                U_in = U[:, 1:n]
                A_sub = U_in' * Dk * U_in
                phi_perm = ryser(A_sub)

                # RHS: exact enumeration over all outputs (incl. collisions)
                phi_exact = exact_char_func(U, n, m, omega, N, k)

                @test abs(phi_perm - phi_exact) < 1e-9
            end
        end
    end

    # Special case: k=0, φ(0) = Per(U_in† U_in) = Per(I_n) = 1
    @testset "k=0: φ(0) = 1 (total probability)" begin
        for (n, m) in [(2, 4), (3, 6), (4, 8)]
            U = RandHaar(m).U
            U_in = U[:, 1:n]
            A = U_in' * U_in
            @test norm(A - I(n)) < 1e-12
            @test abs(ryser(A) - 1.0) < 1e-10
        end
    end
end


# ═══════════════════════════════════════════════════════════════════════
# TEST 6: CDF via DFT inversion (exact, no Monte Carlo)
# S(x0) = (1/N) sum_{k=0}^{N-1} φ(k) G_N(k, x0) (Eq. 11)
# ═══════════════════════════════════════════════════════════════════════
@testset "CDF via exact DFT inversion" begin

    for (n, m) in [(1, 3), (2, 4), (2, 5), (3, 4)]
        @testset "n=$n, m=$m" begin
            Random.seed!(200 + n + 10m)
            U = RandHaar(m).U
            base = n + 1
            omega = [Float64(base^(i - 1)) for i in 1:m]
            N = n * base^(m - 1) + 1

            all_out = enumerate_all_outputs(m, n)
            f_vals = [Int(sum(omega[i] * s[i] for i in 1:m)) for s in all_out]
            test_x0s = sort(unique(f_vals))

            for x0 in test_x0s
                # Exact CDF
                S_exact = exact_cdf(U, n, m, omega, x0)

                # DFT inversion: (1/N) sum_k φ(k) G_N(k, x0)
                S_dft = 0.0
                for k in 0:N-1
                    Dk = diag_Dk(k, omega, N)
                    U_in = U[:, 1:n]
                    phi_k = ryser(U_in' * Dk * U_in)
                    S_dft += real(phi_k * G_N(k, x0, N))
                end
                S_dft /= N

                @test abs(S_dft - S_exact) < 1e-8
            end
        end
    end
end


# ═══════════════════════════════════════════════════════════════════════
# TEST 7: Estimator range bound |Z(k, x0)| = O(log N) (Proposition 3)
# ═══════════════════════════════════════════════════════════════════════
@testset "Estimator range bound (Proposition 3)" begin

    Random.seed!(42)
    m = 6
    n = 2
    U = RandHaar(m).U
    base = n + 1
    omega = [Float64(base^(i - 1)) for i in 1:m]
    N = n * base^(m - 1) + 1

    M = 10_000
    max_Z = 0.0
    for _ in 1:M
        z = Z_sample(U, omega, N, N ÷ 2, n)
        max_Z = max(max_Z, abs(z))
    end

    bound = 10 * log(N)
    @test max_Z < bound
end


# ═══════════════════════════════════════════════════════════════════════
# TEST 8: Monte Carlo CDF estimation accuracy
# ═══════════════════════════════════════════════════════════════════════
@testset "Monte Carlo CDF estimation" begin

    # 8a. Exact comparison for n=2, m=4
    @testset "Exact CDF comparison (n=2, m=4)" begin
        Random.seed!(42)
        n, m = 2, 4
        U = RandHaar(m).U
        base = n + 1
        omega = [Float64(base^(i - 1)) for i in 1:m]
        N = n * base^(m - 1) + 1
        M = 100_000

        all_out = enumerate_all_outputs(m, n)
        f_vals = [Int(sum(omega[i] * s[i] for i in 1:m)) for s in all_out]
        test_x0s = sort(unique(f_vals))

        for x0 in test_x0s
            S_exact = exact_cdf(U, n, m, omega, x0)
            S_mc = estimate_S(U, omega, N, x0, n, M)
            @test abs(S_mc - S_exact) < 0.05
        end
    end

    # 8b. Exact comparison for n=1, m=4
    @testset "Exact CDF comparison (n=1, m=4)" begin
        Random.seed!(123)
        n, m = 1, 4
        U = RandHaar(m).U
        base = n + 1
        omega = [Float64(base^(i - 1)) for i in 1:m]
        N = n * base^(m - 1) + 1
        M = 100_000

        for x0 in [1, 3, 7, N - 1]
            S_exact = exact_cdf(U, n, m, omega, x0)
            S_mc = estimate_S(U, omega, N, x0, n, M)
            @test abs(S_mc - S_exact) < 0.05
        end
    end

    # 8c. Exact comparison for n=3, m=5
    @testset "Exact CDF comparison (n=3, m=5)" begin
        Random.seed!(777)
        n, m = 3, 5
        U = RandHaar(m).U
        base = n + 1
        omega = [Float64(base^(i - 1)) for i in 1:m]
        N = n * base^(m - 1) + 1
        M = 100_000

        all_out = enumerate_all_outputs(m, n)
        f_vals = sort(unique([Int(sum(omega[i] * s[i] for i in 1:m)) for s in all_out]))
        test_x0s = f_vals[round.(Int, range(1, length(f_vals), length=5))]

        for x0 in test_x0s
            S_exact = exact_cdf(U, n, m, omega, x0)
            S_mc = estimate_S(U, omega, N, x0, n, M)
            @test abs(S_mc - S_exact) < 0.05
        end
    end

    # 8d. Monotonicity: S(x0) ≤ S(x0 + 1) (statistically, up to noise)
    @testset "CDF is approximately non-decreasing" begin
        Random.seed!(42)
        n, m = 2, 5
        U = RandHaar(m).U
        base = n + 1
        omega = [Float64(base^(i - 1)) for i in 1:m]
        N = n * base^(m - 1) + 1
        M = 50_000

        prev = -Inf
        for x0 in [0, 5, 10, 15, 20, N - 1]
            S = estimate_S(U, omega, N, x0, n, M)
            @test S > prev - 0.05
            prev = S
        end
    end

    # 8e. S(N-1) ≈ 1 (total probability)
    @testset "S(N-1) ≈ 1" begin
        Random.seed!(42)
        n, m = 2, 4
        U = RandHaar(m).U
        base = n + 1
        omega = [Float64(base^(i - 1)) for i in 1:m]
        N = n * base^(m - 1) + 1
        M = 100_000

        S_total = estimate_S(U, omega, N, N - 1, n, M)
        @test abs(S_total - 1.0) < 0.05
    end
end


# ═══════════════════════════════════════════════════════════════════════
# TEST 9: Binning and unranking
# ═══════════════════════════════════════════════════════════════════════
@testset "Binning and unranking" begin

    # 9a. bin_edges partitions correctly
    @testset "Bins partition correctly" begin
        for (S, d) in [(10, 3), (84, 4), (100, 7), (56, 5), (binomial(8, 3), 4)]
            edges = bin_edges(S, d)
            @test edges[1] == 0
            @test edges[end] == S
            for j in 1:d
                @test edges[j + 1] > edges[j]
            end
        end
    end

    # 9b. Bin widths are ⌊S/d⌋ or ⌊S/d⌋+1
    @testset "Bin width correctness" begin
        for (S, d) in [(15, 3), (84, 4), (252, 6), (100, 7)]
            edges = bin_edges(S, d)
            base_w = S ÷ d
            for j in 1:d
                w = edges[j + 1] - edges[j]
                @test w == base_w || w == base_w + 1
            end
        end
    end

    # 9c. Non-collision state enumeration (legacy)
    @testset "Non-collision state count" begin
        for (m, n) in [(4, 2), (6, 3), (8, 2), (9, 3), (5, 4)]
            states = _enumerate_noncollision_states(m, n)
            @test length(states) == binomial(m, n)
            for s in states
                @test all(x -> x ∈ [0, 1], s)
                @test sum(s) == n
                @test length(s) == m
            end
        end
    end

    # 9d. unrank_composition: all compositions have correct sum
    @testset "unrank_composition: correct sum" begin
        for (n, m) in [(2, 3), (2, 4), (3, 4), (3, 5), (1, 5)]
            total = binomial(n + m - 1, n)
            for rank in 0:total-1
                s = unrank_composition(rank, n, m)
                @test sum(s) == n
                @test length(s) == m
                @test all(x -> x >= 0, s)
            end
        end
    end

    # 9e. unrank_composition: f-values are strictly increasing
    @testset "unrank_composition: f-values strictly increasing" begin
        for (n, m) in [(2, 3), (2, 4), (3, 4), (1, 5), (3, 5)]
            base = n + 1
            total = binomial(n + m - 1, n)
            prev_f = -1
            for rank in 0:total-1
                s = unrank_composition(rank, n, m)
                fv = f_value(s, base)
                @test fv > prev_f
                prev_f = fv
            end
        end
    end

    # 9f. unrank_composition: covers all states (matches enumeration)
    @testset "unrank_composition: covers all states" begin
        for (n, m) in [(2, 3), (2, 4), (3, 4)]
            base = n + 1
            total = binomial(n + m - 1, n)

            # Enumerate all outputs and sort by f-value
            all_out = enumerate_all_outputs(m, n)
            sort!(all_out, by=s -> f_value(s, base))

            @test length(all_out) == total
            for rank in 0:total-1
                s = unrank_composition(rank, n, m)
                @test s == all_out[rank + 1]
            end
        end
    end

    # 9g. f_value: matches manual computation
    @testset "f_value correctness" begin
        @test f_value([2, 0, 0], 3) == 2          # 2*1 + 0*3 + 0*9
        @test f_value([1, 1, 0], 3) == 4           # 1*1 + 1*3 + 0*9
        @test f_value([0, 2, 0], 3) == 6           # 0*1 + 2*3 + 0*9
        @test f_value([0, 0, 2], 3) == 18          # 0*1 + 0*3 + 2*9
        @test f_value([1, 0, 1], 3) == 10          # 1*1 + 0*3 + 1*9
        @test f_value([0, 0, 0, 5], 6) == 5 * 216  # 5 * 6^3
    end

    # 9h. Collision-full f-values are all distinct (injectivity of base-(n+1) encoding)
    @testset "f-values of all states are distinct" begin
        for (n, m) in [(2, 3), (2, 4), (3, 4), (3, 5)]
            base = n + 1
            all_out = enumerate_all_outputs(m, n)
            f_vals = [f_value(s, base) for s in all_out]
            @test length(unique(f_vals)) == length(f_vals)
        end
    end
end


# ═══════════════════════════════════════════════════════════════════════
# TEST 10: Optimized vs legacy path consistency
# ═══════════════════════════════════════════════════════════════════════
@testset "Optimized vs legacy consistency" begin

    # The fused Glynn estimator should produce the same result as the
    # legacy glynn_single on the explicitly formed matrix.
    @testset "Fused Glynn matches legacy" begin
        Random.seed!(42)
        for _ in 1:50
            n = rand(2:5)
            m = n + rand(1:4)
            U = RandHaar(m).U
            U_in = Matrix{ComplexF64}(U[:, 1:n])
            base = n + 1
            N = n * base^(m - 1) + 1
            k = rand(0:min(N - 1, 100))
            x = rand([-1, 1], n)

            # Legacy: form A_sub explicitly
            omega = [Float64(base^(i - 1)) for i in 1:m]
            Dk = diag_Dk(k, omega, N)
            A_sub = U_in' * Dk * U_in
            gly_legacy = glynn_single(A_sub, x)

            # Fused
            buf = Vector{ComplexF64}(undef, m)
            gly_fused = _glynn_fused!(buf, U_in, k, base, N, N ÷ 2, x)

            @test abs(gly_fused - gly_legacy) < 1e-10 * max(1.0, abs(gly_legacy))
        end
    end
end


# ═══════════════════════════════════════════════════════════════════════
# TEST 11: Full pipeline — bin probabilities vs exact (small instance)
# Now using collision-full regime with (n+1)^k weights.
# ═══════════════════════════════════════════════════════════════════════
@testset "Full pipeline: bin probabilities" begin

    # n=2, m=4, d=3 — small enough for exact computation
    @testset "Exact bin probability comparison (n=2, m=4)" begin
        Random.seed!(42)
        n, m, d = 2, 4, 3
        U = RandHaar(m).U            # package convention U[input, output]
        # The local reference helpers (exact_cdf / exact_probability) are written in
        # the standard A[output, input] orientation, whereas find_most_probable_bin
        # takes the package convention and transposes internally — so the reference
        # must be fed the transpose of the same matrix.
        Uout = permutedims(U)
        base = n + 1
        omega = [Float64(base^(i - 1)) for i in 1:m]
        N = n * base^(m - 1) + 1
        total_states = binomial(n + m - 1, n)

        # All states sorted by f-value (via unranking)
        edges = bin_edges(total_states, d)

        # Exact CDF at bin edges
        exact_bin_probs = Float64[]
        S_prev = 0.0
        for j in 1:d
            if edges[j + 1] == total_states
                S_curr = 1.0
            else
                s = unrank_composition(edges[j + 1] - 1, n, m)
                x_j = f_value(s, base)
                S_curr = exact_cdf(Uout, n, m, omega, x_j)
            end
            push!(exact_bin_probs, S_curr - S_prev)
            S_prev = S_curr
        end

        # Monte Carlo
        j_star_mc, mc_bin_probs = find_most_probable_bin(U, m, n, d; M=100_000)
        j_star_exact = argmax(exact_bin_probs)

        for j in 1:d
            @test abs(mc_bin_probs[j] - exact_bin_probs[j]) < 0.05
        end
        # The argmax is only a meaningful assertion when the top two exact bins are
        # separated by more than the MC tolerance above.  estimate_S draws inside a
        # Threads.@threads loop with task-local RNGs, so Random.seed! does not pin
        # these samples and the result varies with JULIA_NUM_THREADS; asserting a
        # hard equality across a near-tie makes the test flaky rather than strict.
        sorted_exact = sort(exact_bin_probs; rev=true)
        if sorted_exact[1] - sorted_exact[2] > 0.10
            @test j_star_mc == j_star_exact
        else
            @test mc_bin_probs[j_star_exact] > maximum(mc_bin_probs) - 0.10
        end
    end

    # Dilute regime: n=3, m=9 — bin probs should sum ≈ 1
    @testset "Dilute regime (n=3, m=9)" begin
        Random.seed!(123)
        n, m, d = 3, 9, 4
        U = RandHaar(m).U

        _, bin_probs = find_most_probable_bin(U, m, n, d; M=20_000)

        @test sum(bin_probs) > 0.85
        @test all(p -> p > -0.05, bin_probs)
    end
end


# ═══════════════════════════════════════════════════════════════════════
# TEST 12: Convergence rate — error ∝ 1/√M
# ═══════════════════════════════════════════════════════════════════════
@testset "Convergence rate" begin

    Random.seed!(42)
    n, m = 2, 4
    U = RandHaar(m).U
    base = n + 1
    omega = [Float64(base^(i - 1)) for i in 1:m]
    N = n * base^(m - 1) + 1
    x0 = 5
    S_exact = exact_cdf(U, n, m, omega, x0)

    Ms = [1_000, 10_000, 100_000]
    errors = Float64[]
    for M_val in Ms
        errs = Float64[]
        for _ in 1:30
            S_mc = estimate_S(U, omega, N, x0, n, M_val)
            push!(errs, abs(S_mc - S_exact))
        end
        push!(errors, mean(errs))
    end

    @test errors[2] < errors[1]
    @test errors[3] < errors[1]
    ratio = errors[1] / errors[3]
    @test ratio > 2.0
end

# ═══════════════════════════════════════════════════════════════════════
# TEST 13: Large-N Int128 / BigInt path
# ═══════════════════════════════════════════════════════════════════════
@testset "Large-N integer types" begin

    # 13a. Int128 path: n=15 overflows Int64 (N ≈ 1.7e19)
    @testset "Int128 path (n=15, m=21)" begin
        Random.seed!(42)
        n, m, d = 15, 21, 4
        U = RandHaar(m).U
        N = _compute_N(n, n + 1, m)
        @test N isa Int128
        @test N > typemax(Int64)

        ctx = SamplingContext(N)
        @test ctx.N isa Int128
        k = sample_fourier_mode(ctx)
        @test k isa Int128

        fv = f_value(unrank_composition(0, n, m), Int128(n + 1))
        @test fv isa Int128

        j_star, bin_probs = find_most_probable_bin(U, m, n, d; M=500)
        @test 1 <= j_star <= d
        @test length(bin_probs) == d
    end

    # 13b. Auto-selection: small (n,m) stays Int64
    @testset "Int64 path for small instances" begin
        N = _compute_N(3, 4, 10)
        @test N isa Int
    end

    # 13c. Fused Glynn consistency with Int128 types
    @testset "Fused Glynn matches legacy (Int128)" begin
        Random.seed!(42)
        n, m = 3, 5
        U = RandHaar(m).U
        U_in = Matrix{ComplexF64}(U[:, 1:n])
        base128 = Int128(n + 1)
        N128 = Int128(n) * base128^(m - 1) + Int128(1)
        k128 = Int128(7)
        x = rand([-1, 1], n)

        # Legacy
        omega = [Float64(base128^(i - 1)) for i in 1:m]
        Dk = diag_Dk(k128, omega, N128)
        gly_legacy = glynn_single(U_in' * Dk * U_in, x)

        # Fused with Int128 types
        buf = Vector{ComplexF64}(undef, m)
        gly_fused = _glynn_fused!(buf, U_in, k128, base128, N128, N128 ÷ 2, x)

        @test abs(gly_fused - gly_legacy) < 1e-10 * max(1.0, abs(gly_legacy))
    end
end


# ═══════════════════════════════════════════════════════════════════════
# TEST 14: Fast Clifford & Clifford sampler + z_samples
# ═══════════════════════════════════════════════════════════════════════
@testset "Clifford & Clifford sampler and z_samples" begin

    # 14a. Every cc_sample! output is a valid collision-full configuration
    @testset "cc_sample! output validity" begin
        Random.seed!(7)
        n, m = 3, 6
        U = RandHaar(m).U
        ws = CCSamplerWorkspace(U, n)
        for _ in 1:2000
            out = cc_sample!(ws)
            @test length(out) == n          # n photons placed
            @test all(md -> 1 <= md <= m, out)   # valid modes
            @test issorted(out)             # returned sorted
        end
    end

    # 14b. Empirical PMF converges to the exact boson-sampling PMF
    #      Pr[s] = |Per(A[s_modes, 1:n])|² / ∏ s_j!  (collisions included), where
    #      A = transpose(U) because CCSamplerWorkspace takes U[input, output].
    @testset "empirical PMF matches |Per|²/∏s_j! (n=2, m=3)" begin
        Random.seed!(42)
        n, m = 2, 3
        U = RandHaar(m).U
        Uout = permutedims(U)                      # reference helper wants [output, input]
        states = enumerate_all_outputs(m, n)       # all C(n+m-1,n) collision-full outputs
        idx = Dict(s => i for (i, s) in enumerate(states))

        K = 200_000
        ws = CCSamplerWorkspace(U, n)
        counts = zeros(Int, length(states))
        occ = zeros(Int, m)
        for _ in 1:K
            out = cc_sample!(ws)
            fill!(occ, 0)
            for md in out; occ[md] += 1; end
            counts[idx[copy(occ)]] += 1
        end
        emp = counts ./ K
        exact = [exact_probability(Uout, n, s) for s in states]

        @test isapprox(sum(exact), 1.0; atol=1e-10)   # exact PMF is normalized
        @test sum(emp) ≈ 1.0                          # sampler always lands on a state
        # total-variation distance ≪ 1: per-bin error ~ 1/√K ≈ 2e-3
        @test 0.5 * sum(abs.(emp .- exact)) < 0.01
    end

    # 14b'. Convention regression: the sampler must agree with the package's own
    #       compute_probability! on the SAME interferometer object.  The local
    #       exact_probability helper above cannot catch a transposed U (it shares
    #       the orientation with the kernel and RandHaar is transpose-invariant in
    #       law); compute_probability! is the independent source of truth, since
    #       scattering_matrix() indexes U[index_input, index_output].
    #       Measured: TVD 0.0016 with the correct convention vs 0.339 transposed.
    @testset "convention matches compute_probability!" begin
        Random.seed!(11)
        n, m = 2, 3
        interf = RandHaar(m)
        states = enumerate_all_outputs(m, n)
        idx = Dict(s => i for (i, s) in enumerate(states))

        pkg = map(states) do s
            ev = Event(Input{Bosonic}(first_modes(n, m)),
                       FockDetection(ModeOccupation(s)), interf)
            compute_probability!(ev)
            ev.proba_params.probability
        end
        @test isapprox(sum(pkg), 1.0; atol=1e-8)

        K = 200_000
        ws = CCSamplerWorkspace(interf, n)   # Interferometer method: no manual transpose
        counts = zeros(Int, length(states))
        occ = zeros(Int, m)
        for _ in 1:K
            out = cc_sample!(ws)
            fill!(occ, 0)
            for md in out; occ[md] += 1; end
            counts[idx[copy(occ)]] += 1
        end
        emp = counts ./ K
        # ≈2e-3 per-bin MC noise; a transposed convention lands two orders up.
        @test 0.5 * sum(abs.(emp .- pkg)) < 0.02
    end

    # 14c. mean(z_samples) reproduces estimate_S — both average the same iid Z
    @testset "mean(z_samples) ≈ estimate_S" begin
        Random.seed!(2024)
        n, m = 3, 5
        U = RandHaar(m).U
        base = n + 1
        N = _compute_N(n, base, m)
        ctx = SamplingContext(N)
        U_in = U[:, 1:n]
        x0 = N ÷ 3

        K = 400_000
        zs = z_samples(U_in, base, N, x0, n, K, ctx)
        @test length(zs) == K
        @test all(isfinite, zs)

        S_direct = estimate_S(U_in, base, N, x0, n, K, ctx)
        # two independent K-sample means of the same estimator: agree to a few SE
        se = std(zs) / sqrt(K)
        @test abs(mean(zs) - S_direct) < 8 * se + 1e-3
        @test -0.2 < mean(zs) < 1.2                    # sane CDF range
    end
end

end  # @testset "one-way function (boson sampling OWF)"
