# Cost of ONE value of the cumulative function S(x0) at n=m=100, with a VERIFIED error.
#
#   Run:  JULIA_NUM_THREADS=24 julia --project=. check_n100_1percent.jl
#
# Cost model:   M_required = Var(Z)/ε²   and   per-sample time (O(m·n), thread-divided).
#
# ε is the 1-σ Monte-Carlo STANDARD ERROR of Ŝ(x0), not a confidence bound: even at the
# target ε, ~32% of estimates deviate by more than ε (a 95% bound needs ~4× the samples).
#
# Two x0 are timed because the per-sample variance Var(Z) — and hence the cost — varies
# strongly across the f-value distribution:
#
#   x0 = N/1000  — the MEDIAN (S ≈ 0.5).  The f-values are concentrated at small f, so
#                  S climbs 0→1 over x0/N ∈ [1e-6, 0.1].  x0 = N/1000 = 0.1·(n+1)^(m-1)
#                  forces the leading digit s_m = 0, hence S(N/1000) = P(last mode empty)
#                  ≈ 0.5.  This is a LOW-variance point — NOT the worst MC case.
#
#   x0 = N/2     — the bulk (S ≈ 1).  Var(Z) is ~4-6× larger here, so this is close to the
#                  worst-case cost of a single CDF value.  The "literal middle" x0 = N/2
#                  has S ≈ 1, NOT 0.5.
#
# Var(Z) calibrated from a few short runs is noisy, so the predicted M_req can run a bit
# optimistic.  The script therefore VERIFIES with a second, independent Var(Z) reading:
# n_verify estimates at M_ver = min(M_req, cap).  Since |Z| = O(log N) is bounded, Var(Z)
# is well estimated without spending the full M_req per run — the achieved SD is rescaled
# to M_req via the exact mean-variance law Var(Ŝ_M) = Var(Z)/M.  This avoids the (dominant)
# n_verify × M_req verification cost, which at the bulk x0 = N/2 alone runs ~6+ min.

using LinearAlgebra, Statistics, Printf
using BosonSampling  # OWF estimator is compiled in & auto-exported

fmt(s) = s < 1 ? @sprintf("%.0f ms", s*1000) : @sprintf("%.2f s", s)

n = m = 100; ε = 0.01
base = n + 1
N = _compute_N(n, base, n); T = typeof(N)
U_in = Matrix{ComplexF64}(RandHaar(m).U[:, 1:n])
ctx = SamplingContext(N)

# Calibrate Var(Z) at x0 = q·N, predict M_req for a 1-σ error ε, then CHECK the prediction
# by running n_verify independent estimates at M_req and reporting their actual SD.
function bench(label, q; n_cal=20, M_cal=200_000, n_verify=12, M_verify_cap=1_000_000)
    x0 = T(round(BigInt, big(N) * q))
    estimate_S(U_in, T(base), N, x0, n, 2000, ctx)                       # warmup / JIT

    # 1) calibrate per-sample time and Var(Z)
    t_cal = @elapsed Svals = [estimate_S(U_in, T(base), N, x0, n, M_cal, ctx) for _ in 1:n_cal]
    t_per_sample = t_cal / (n_cal * M_cal)
    VarZ  = var(Svals) * M_cal
    M_req = ceil(Int, VarZ / ε^2)

    # 2) empirical verification: a second, independent Var(Z) reading from n_verify estimates
    #    at M_ver = min(M_req, cap).  Var(Ŝ_M) = Var(Z)/M is exact, so we rescale the spread
    #    to the SD M_req actually achieves — no need to pay n_verify × M_req.
    M_ver    = min(M_req, M_verify_cap)
    Sver     = [estimate_S(U_in, T(base), N, x0, n, M_ver, ctx) for _ in 1:n_verify]
    VarZ_ver = var(Sver) * M_ver
    sd_hit   = sqrt(VarZ_ver / M_req)            # SD that the predicted M_req achieves
    M_true   = ceil(Int, VarZ_ver / ε^2)         # M to actually reach ε, from verified Var(Z)

    @printf("\n[%s]   x0/N ≈ %.0e\n", label, q)
    @printf("  S(x0) ≈ %.3f    Var(Z) ≈ %.0f    t/sample = %.2f µs (wall, %d threads)\n",
            mean(Svals), VarZ, t_per_sample * 1e6, Threads.nthreads())
    @printf("  predicted: M_req = %d  ⇒  one S(x0) ≈ %s\n", M_req, fmt(M_req * t_per_sample))
    @printf("  VERIFIED (%d runs @ M_ver=%d): Var(Z) ≈ %.0f  ⇒  SD@M_req = %.4f  (target ε = %.4f)  →  %.2f%%\n",
            n_verify, M_ver, VarZ_ver, sd_hit, ε, 100 * sd_hit)
    @printf("  to actually hit ε=%.0f%%: M ≈ %d  ⇒  one S(x0) ≈ %s\n",
            ε * 100, M_true, fmt(M_true * t_per_sample))
end

@printf("threads = %d   (N has %d digits)\n", Threads.nthreads(), ndigits(N))
bench("median  S≈0.5  (low variance)", 1e-3)
bench("bulk    S≈1    (~worst case)",  0.5)
