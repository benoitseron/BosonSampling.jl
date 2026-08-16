# ════════════════════════════════════════════════════════════════════
# Validation: does the Monte-Carlo estimator reproduce the EXACT numbers?
# ════════════════════════════════════════════════════════════════════
#
# Panel (d) of the paper figure turns a per-shot variance σ² into a sample
# budget M = (z/ε)²·σ². Two things must be true for that to be meaningful:
#   1. the MC-Fourier estimator Ŝ(x₀) is UNBIASED for the true S(x₀);
#   2. spending the prescribed M samples really does land you within ±ε
#      of the truth ~95% of the time.
# Here we check both against ground truth obtained by DIRECT SUMMATION of
# every output probability — feasible at a small n where the whole output
# space can be enumerated.
#
#   (a) S(x₀) across thresholds : sampled mean ± 95% CI  vs  exact CDF.
#   (b) budget / coverage       : with M = (z/ε)²σ² samples, what fraction
#                                 of estimates fall within ±ε of exact S?

using LinearAlgebra, Random, Printf, Statistics
using BosonSampling  # OWF estimator is compiled in & auto-exported
using Permanents: ryser
using Combinatorics: multiexponents
using Plots, LaTeXStrings

default(
    size          = (720, 460),
    bottom_margin = 6 * Plots.mm, left_margin = 8 * Plots.mm,
    right_margin  = 4 * Plots.mm, top_margin  = 4 * Plots.mm,
    titlefont = font("sans-serif", 11), guidefont = font("sans-serif", 10),
    tickfont  = font("sans-serif",  9), legendfont = font("sans-serif", 8),
    framestyle = :box, grid = true, gridalpha = 0.25,
)
COLOR_MCF = RGB(0.85, 0.40, 0.10)   # Monte-Carlo estimates
COLOR_EX  = RGB(0.20, 0.45, 0.70)   # exact ground truth

# ── exact ground truth helpers (lifted from paper_figure.jl / runtests.jl) ──
# Weak compositions of n photons into m modes = Combinatorics.multiexponents.
enumerate_all_outputs(m::Int, n::Int) = collect(multiexponents(m, n))
# Pr[s] = |Per(U[out_cols, 1:n])|² / ∏ sⱼ!  for input |1ⁿ0^{m-n}⟩
function exact_probability(U::AbstractMatrix, n::Int, output::Vector{Int})
    out_cols = Int[]
    for j in eachindex(output), _ in 1:output[j]; push!(out_cols, j); end
    return abs2(ryser(U[out_cols, 1:n])) / prod(factorial(s) for s in output)
end

# ── setup: one fixed device, small enough to enumerate exactly ──
n    = 6                              # exact enumeration: C(n+m-1, n) output states
m    = round(Int, 2.1n)               # same shallow regime as panel (d): m = 2.1n
base = n + 1
N    = _compute_N(n, base, m); B = typeof(N)(base)
ctx  = SamplingContext(N)
ε, z = 0.05, 1.96
Random.seed!(20260713)
V    = RandHaar(m).U                  # a fixed Haar interferometer (the "device")
U_in = V[:, 1:n]

# ── EXACT S(x₀) by direct summation of all output probabilities ──
states = enumerate_all_outputs(m, n)
fvals  = [f_value(s, B) for s in states]
probs  = [exact_probability(V, n, s) for s in states]
@assert isapprox(sum(probs), 1.0; atol = 1e-9) "exact probabilities must sum to 1"
perm    = sortperm(fvals)
fsorted = fvals[perm]                 # output f-values in increasing order
cdf     = cumsum(probs[perm])         # exact S at each rank boundary: S = P(f ≤ fsorted[r])
total   = length(states)
@printf("n=%d  m=%d  →  %d output states enumerated exactly\n\n", n, m, total)

# ════════════════════════════════════════════════════════════════════
# (a) sampled S(x₀) vs exact, across thresholds — BOTH estimators
# ════════════════════════════════════════════════════════════════════
K   = 20000                           # Z-shots per threshold (MC-Fourier)
Kcc = 20000                           # Clifford-Clifford shots (drawn once, reused)

# Direct quantum sampling is threshold-independent: draw Kcc shots ONCE, encode each
# to its base-B f-value, then every S_cc(x₀) is just a fraction below x₀.
# This script works in the kernel orientation V[output, input] (see exact_probability
# above); CCSamplerWorkspace is user-facing and takes the package convention
# V[input, output], so hand it the transpose to stay on the same device.
ws  = CCSamplerWorkspace(permutedims(V), n); pw = [B^(k - 1) for k in 1:m]
fcc = [sum(pw[md] for md in cc_sample!(ws)) for _ in 1:Kcc]

ranks_a = unique(round.(Int, range(1, total, length = 13)))
qa, Sex = Float64[], Float64[]
Smc, Sci = Float64[], Float64[]       # MC-Fourier mean ± CI
Scc, Scci = Float64[], Float64[]      # direct CC fraction ± CI
println("─── (a) S(x₀): sampled vs exact ───")
println("  quantile  S_exact   S_MC±CI   pull   S_CC±CI   pull")
for r in ranks_a
    x0 = fsorted[r]
    Z  = z_samples(U_in, B, N, x0, n, K, ctx)     # K raw MC-Fourier single-shot estimates
    scc = count(≤(x0), fcc) / Kcc                 # Bernoulli(S) fraction below x₀
    push!(qa, r / total); push!(Sex, cdf[r])
    push!(Smc, mean(Z));  push!(Sci, z * std(Z) / sqrt(K))
    push!(Scc, scc);      push!(Scci, z * sqrt(scc * (1 - scc) / Kcc))
    # pull = (estimate − exact)/CI; undefined for CC when its CI collapses to 0
    # (scc ∈ {0,1} at the extreme thresholds — a Bernoulli with no spread).
    ccpull = Scci[end] > 1e-6 ? (Scc[end] - cdf[r]) / Scci[end] : 0.0
    @printf("   %.3f    %.4f   %.4f   %+.2f   %.4f   %+.2f\n", r / total, cdf[r],
        Smc[end], (Smc[end] - cdf[r]) / Sci[end], Scc[end], ccpull)
end

p_a = plot((1:total) ./ total, cdf;                                   # exact CDF (line)
    color = COLOR_EX, lw = 2, label = "exact (direct summation)",
    xlabel = "threshold quantile " * L"r / N_\mathrm{states}",
    ylabel = L"S(x_0) = P(f \leq x_0)", title = "(a) both estimators vs exact CDF",
    legend = :topleft, dpi = 800)
scatter!(p_a, qa, Smc; yerror = Sci, color = COLOR_MCF, marker = :circle,
    markersize = 5, label = "MC-Fourier (classical)")
scatter!(p_a, qa, Scc; yerror = Scci, color = RGB(0.20, 0.60, 0.35), marker = :diamond,
    markersize = 5, label = "direct sampling (quantum)")

# ════════════════════════════════════════════════════════════════════
# (b) budget / coverage at the median threshold (worst case, S ≈ 0.5)
# ════════════════════════════════════════════════════════════════════
# Spend exactly the panel-(d) budget M = (z/ε)²·σ² and ask: do the estimates
# actually fall within ±ε of the exact S the promised 95% of the time?
r_med = cld(total, 2)
x0m, S0 = fsorted[r_med], cdf[r_med]
σ2    = var(z_samples(U_in, B, N, x0m, n, 200_000, ctx))  # reference per-shot variance
Mstar = ceil(Int, (z / ε)^2 * σ2)                         # prescribed sample budget
T     = 1000                                              # independent trials
resid = [estimate_S(U_in, B, N, x0m, n, Mstar, ctx) - S0 for _ in 1:T]
cover = count(r -> abs(r) ≤ ε, resid) / T                # empirical 95%-CI coverage

println("\n─── (b) budget / coverage at median threshold ───")
@printf("  exact S(x₀)      = %.4f\n", S0)
@printf("  per-shot σ²       = %.3f  →  M = (z/ε)²σ² = %d samples\n", σ2, Mstar)
@printf("  empirical coverage of ±ε=%.2f : %.1f%%   (target 95%%)\n", ε, 100cover)

p_b = histogram(resid; bins = 40, normalize = :pdf, color = COLOR_MCF,
    alpha = 0.55, label = "estimates (M = $Mstar)",
    xlabel = L"\hat{S} - S_\mathrm{exact}", ylabel = "density",
    title = "(b) coverage: " * @sprintf("%.0f%% within ±ε", 100cover),
    legend = :topright, dpi = 800)
σŜ = ε / z                                                # predicted std of Ŝ at M = Mstar
xs = range(minimum(resid), maximum(resid), length = 200)
plot!(p_b, xs, exp.(-(xs .^ 2) ./ (2σŜ^2)) ./ (σŜ * sqrt(2π));
    color = COLOR_EX, lw = 2, label = "predicted " * L"\mathcal{N}(0,(ε/z)^2)")
vline!(p_b, [-ε, ε]; color = :black, ls = :dash, lw = 1, label = L"\pm ε")

fig = plot(p_a, p_b; layout = (1, 2), size = (1200, 460),
    bottom_margin = 7 * Plots.mm, left_margin = 9 * Plots.mm)
savefig(fig, joinpath(@__DIR__, "validate_sampling.pdf"))
savefig(fig, joinpath(@__DIR__, "validate_sampling.png"))
println("\nsaved validate_sampling.pdf and validate_sampling.png")
