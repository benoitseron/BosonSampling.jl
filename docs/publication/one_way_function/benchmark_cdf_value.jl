# Benchmark: cost of generating ONE value of the cumulative function S(x0).
#
# S(x0) is computed by estimate_S (src/main.jl:394) via M Monte-Carlo samples,
# each an O(m·n) fused-Glynn evaluation. This script isolates the wall-clock
# cost of a single estimate_S call and confirms the two scaling laws:
#   (1) linear in M           (the MC budget)
#   (2) linear in m·n         (per-sample, no permanent)
#
# Run with:  julia --project=. benchmark_cdf_value.jl
# Threads:   JULIA_NUM_THREADS=N julia --project=. benchmark_cdf_value.jl

using LinearAlgebra, Random, Printf, Statistics
using BosonSampling

const NT = Threads.nthreads()
fmt(s) = s < 1e-3 ? @sprintf("%.0f µs", s*1e6) :
         s < 1    ? @sprintf("%.2f ms", s*1e3) :
                    @sprintf("%.3f s", s)

# Build a problem instance and return the args for one estimate_S call.
function setup(n, m)
    base = n + 1
    N = _compute_N(n, base, m)
    T = typeof(N)
    U = RandHaar(m).U
    U_in = Matrix{ComplexF64}(U[:, 1:n])
    ctx = SamplingContext(N)
    x0 = N ÷ 2
    return (U_in, T(base), N, T(x0), n, ctx)
end

# Time a single estimate_S call (one CDF value), best of `reps` runs.
function time_one_value(n, m, M; reps=5)
    U_in, base, N, x0, nn, ctx = setup(n, m)
    estimate_S(U_in, base, N, x0, nn, min(M, 2000), ctx)   # warmup / compile
    best = Inf
    for _ in 1:reps
        t = @elapsed estimate_S(U_in, base, N, x0, nn, M, ctx)
        best = min(best, t)
    end
    return best, typeof(N)
end

println("="^64)
println("ONE value of the cumulative function S(x0)  —  estimate_S")
println("Julia threads = $NT   (machine cores: ", Sys.CPU_THREADS, ")")
println("="^64)

# ─── 1. Headline: time for one CDF value at representative sizes ──────
println("\n[1] Wall-clock for ONE S(x0) value  (M = 50,000 samples)")
println("    n    m       N-type        time      µs/sample")
M0 = 50_000
for (n, m) in [(5,15), (10,15), (14,15), (20,21), (25,26), (30,31), (40,40)]
    t, Nt = time_one_value(n, m, M0)
    @printf("   %2d   %2d   %10s   %9s   %7.3f\n", n, m, Nt, fmt(t), t/M0*1e6)
end

# ─── 2. Linear in M ──────────────────────────────────────────────────
println("\n[2] Scaling in M  (fixed n=12, m=15)")
println("        M        time      µs/sample")
n2, m2 = 12, 15
for M in [5_000, 10_000, 20_000, 50_000, 100_000, 200_000]
    t, _ = time_one_value(n2, m2, M)
    @printf("   %8d   %9s   %7.3f\n", M, fmt(t), t/M*1e6)
end

# ─── 3. Linear in m·n (per-sample cost) ──────────────────────────────
println("\n[3] Per-sample cost vs m·n  (M = 30,000)")
println("    n    m    m·n      µs/sample")
M3 = 30_000
cases3 = [(5,8),(5,16),(5,24),(8,15),(12,15),(14,15),(10,11),(13,14),(20,21),(25,26)]
mn_v = Float64[]; us_v = Float64[]
for (n, m) in cases3
    Float64(n)*Float64(n+1)^(m-1) > 1.7e38 && continue   # keep off the BigInt path for clean fit
    t, _ = time_one_value(n, m, M3)
    us = t/M3*1e6
    push!(mn_v, n*m); push!(us_v, us)
    @printf("   %2d   %2d   %4d      %7.3f\n", n, m, n*m, us)
end
# least-squares fit  us ≈ a·(mn) + b
k = length(mn_v)
a = (k*sum(mn_v.*us_v) - sum(mn_v)*sum(us_v)) / (k*sum(mn_v.^2) - sum(mn_v)^2)
b = mean(us_v) - a*mean(mn_v)
@printf("\n   fit:  µs/sample ≈ %.5f·(m·n) + %.3f\n", a, b)
@printf("   => one S(x0) value ≈ M · (%.2e·mn + %.2e) seconds  (÷ threads)\n", a*1e-6, b*1e-6)
