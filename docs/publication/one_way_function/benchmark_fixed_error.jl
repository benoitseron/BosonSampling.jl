# How long to get ONE value of S(x0) at FIXED statistical error ε?
#
# S(x0) is the mean of M i.i.d. samples Z (estimate_S, src/main.jl:394).
# MC standard error = sqrt(Var Z / M).  To reach absolute error ε:
#
#       M_required = Var(Z) / ε²
#       time(one value) = M_required · (per-sample time) / threads
#
# The per-sample TIME is polynomial O(m·n); the open question is whether
# Var(Z) stays manageable as n grows.  We measure it directly at n=m=100.
#
# Run:  JULIA_NUM_THREADS=N julia --project=. benchmark_fixed_error.jl

using LinearAlgebra, Random, Printf, Statistics
using BosonSampling

const NT = Threads.nthreads()
fmt(s) = s < 1e-3 ? @sprintf("%.0f µs", s*1e6) :
         s < 1     ? @sprintf("%.2f ms", s*1e3) :
         s < 60    ? @sprintf("%.2f s", s) :
         s < 3600  ? @sprintf("%.1f min", s/60) :
         s < 86400 ? @sprintf("%.1f h", s/3600) :
                     @sprintf("%.1f days", s/86400)

function setup(n, m)
    base = n + 1
    N = _compute_N(n, base, m)
    T = typeof(N)
    U = RandHaar(m).U
    U_in = Matrix{ComplexF64}(U[:, 1:n])
    ctx = SamplingContext(N)
    return (U_in, T(base), N, n, ctx)
end

# Estimate Var(Z): run `trials` independent estimate_S calls, each averaging
# M_var samples; the spread of the returned means gives Var(Z) = trials-std² · M_var.
function measure(n, m, x0frac; M_time=50_000, M_var=200_000, trials=30)
    U_in, base, N, nn, ctx = setup(n, m)
    T = typeof(N)
    x0 = T(round(BigInt, BigInt(N) * x0frac))

    # per-sample wall time
    estimate_S(U_in, base, N, x0, nn, 2000, ctx)            # warmup
    t = @elapsed estimate_S(U_in, base, N, x0, nn, M_time, ctx)
    t_per_sample = t / M_time                               # already ÷ threads (estimate_S is threaded)

    # variance of the estimator
    Svals = Float64[]
    for _ in 1:trials
        push!(Svals, estimate_S(U_in, base, N, x0, nn, M_var, ctx))
    end
    S̄ = mean(Svals)
    VarZ = var(Svals) * M_var          # Var of one sample
    return (; S̄, VarZ, t_per_sample, N)
end

ε = 0.01
println("="^70)
println("FIXED error ε = $ε (absolute, on S) — required M and wall time")
println("Julia threads = $NT")
println("="^70)

for (n, m) in [(40,40), (100,100)]
    println("\n── n = $n, m = $m ─────────────────────────────")
    println("  x0/N    S(x0)     Var(Z)     M_req=Var/ε²    t/sample      time(1 value)")
    for x0frac in (0.25, 0.5, 0.75)
        r = measure(n, m, x0frac)
        M_req = r.VarZ / ε^2
        time_val = M_req * r.t_per_sample
        @printf("  %.2f   %6.3f   %8.3f   %12s   %9s   %s\n",
                x0frac, r.S̄, r.VarZ, @sprintf("%.2e", M_req),
                fmt(r.t_per_sample), fmt(time_val))
    end
end
