"""
Distribution-correctness tests for the Householder sampler.

The existing `test_householder_integration.jl` and
`test_one_parameter_householder.jl` tests only check photon-number
conservation and vector length. Those pass even when the sampler
produces a wrong distribution. The tests below compare the empirical
sampler distribution against `full_distribution` / `compute_probability!`
and flag the regression that slips past photon conservation.

Run with: `julia --project=. test/test_householder_distribution.jl`
from the package root. Requires `using Test, BosonSampling, Random,
LinearAlgebra, StatsBase`.

KNOWN FAILURE (reported by Ulysse, 2026-04-20): every test below fails
with TVD ≈ 0.5 against the exact distribution, even in the Bosonic
limit `S = ones(n, n)` where the Householder sampler must reduce to
ordinary Clifford. The photon-conservation tests still pass.
"""

using Test
using BosonSampling
using Random
using LinearAlgebra
using StatsBase: Weights

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

"""Histogram of mode-occupation samples into a Dict."""
function histogram(samples::AbstractVector)
    h = Dict{Vector{Int}, Int}()
    for s in samples
        h[s] = get(h, s, 0) + 1
    end
    return h
end

"""TVD between empirical samples and exact probabilities over the
outcome list `outcomes` (each a mode-occupation vector)."""
function empirical_tvd(samples::AbstractVector,
                       outcomes::AbstractVector,
                       probs_exact::AbstractVector)
    N = length(samples)
    h = histogram(samples)
    t = 0.0
    for (k, o) in enumerate(outcomes)
        t += abs(get(h, o, 0) / N - probs_exact[k])
    end
    return t / 2
end

"""Draw N samples by calling `householder_sampler` directly (bypasses
sample!/sample_multiple dispatch), so we know we are hitting the
Householder code path and nothing else."""
function draw_householder(input, interf, N::Integer)
    out = Vector{Vector{Int}}(undef, N)
    for i in 1:N
        ev = Event(input, FockSample(), interf)
        out[i] = householder_sampler(ev).state
    end
    return out
end

"""Exact distribution via full_distribution, returned as
(outcomes, probs) with probs normalised."""
function exact_distribution(input, interf)
    dist = full_distribution(input, interf)
    outcomes = [copy(c.state) for c in dist.counts]
    probs = Float64.(real.(dist.proba))
    probs ./= sum(probs)
    return outcomes, probs
end

# Fixed seed so every run of this file is bit-identical.
SEED = 20260420

# All diagnostics inside one outer testset so that a failure in T2 does
# NOT halt the rest — we want the full picture in every run.
@testset "Householder sampler: distribution correctness" begin

# ---------------------------------------------------------------------------
# T0 — sanity: full_distribution(UserDefinedGramMatrix, S=ones) equals
# full_distribution(Bosonic). This is NOT a sampler test — it verifies
# that the reference distribution we are comparing against is the right
# one. If T0 fails, the bug is in `process_probability_partial`, not the
# sampler.
# ---------------------------------------------------------------------------
@testset "T0: full_distribution(G=ones) == full_distribution(Bosonic)" begin
    Random.seed!(SEED)
    n, m = 3, 5
    interf = RandHaar(m)

    input_B = Input{Bosonic}(first_modes(n, m))
    input_G = Input{UserDefinedGramMatrix}(first_modes(n, m),
                                           ComplexF64.(ones(n, n)))

    _, p_B = exact_distribution(input_B, interf)
    _, p_G = exact_distribution(input_G, interf)
    @test p_B ≈ p_G atol=1e-10      # must be *identical* — same physical state
end

# ---------------------------------------------------------------------------
# T1 — sanity: Clifford on Input{Bosonic} matches the exact distribution.
# This is a passing baseline. If this FAILS, the Clifford sampler itself
# is broken, not the Householder layer — so debug Clifford first.
# ---------------------------------------------------------------------------
@testset "T1 [BASELINE, should pass]: cliffords_sampler(Bosonic) vs exact" begin
    Random.seed!(SEED)
    n, m, N = 3, 5, 8_000
    interf = RandHaar(m)
    input  = Input{Bosonic}(first_modes(n, m))

    outcomes, probs = exact_distribution(input, interf)
    samples = [cliffords_sampler(input=input, interf=interf) for _ in 1:N]
    tvd = empirical_tvd(samples, outcomes, probs)
    println("T1 TVD (Clifford Bosonic) = ", round(tvd, digits=4))
    @test tvd < 0.05    # ≳ statistical noise at N=8k, ≪ ~0.5 failure signature
end

# ---------------------------------------------------------------------------
# T2 — the key failing case. Householder(S=ones) must reduce to Bosonic
# Clifford by construction: S = ones is the Gram matrix of fully
# indistinguishable photons.
# ---------------------------------------------------------------------------
@testset "T2 [FAILS]: householder_sampler(G=ones) vs exact Bosonic" begin
    Random.seed!(SEED)
    n, m, N = 3, 5, 8_000
    interf = RandHaar(m)
    input  = Input{UserDefinedGramMatrix}(first_modes(n, m),
                                          ComplexF64.(ones(n, n)))

    outcomes, probs = exact_distribution(input, interf)
    samples = draw_householder(input, interf, N)
    tvd = empirical_tvd(samples, outcomes, probs)
    println("T2 TVD (Householder G=ones) = ", round(tvd, digits=4))
    @test tvd < 0.05
end

# ---------------------------------------------------------------------------
# T3 — fully-distinguishable limit. S = I. Must reduce to the classical
# sampler's distribution, which is `full_distribution(Distinguishable)`.
# ---------------------------------------------------------------------------
@testset "T3 [FAILS]: householder_sampler(G=I) vs exact Distinguishable" begin
    Random.seed!(SEED)
    n, m, N = 3, 5, 8_000
    interf = RandHaar(m)
    input  = Input{UserDefinedGramMatrix}(first_modes(n, m),
                                          ComplexF64.(Matrix(I, n, n)))

    outcomes, probs = exact_distribution(input, interf)
    samples = draw_householder(input, interf, N)
    tvd = empirical_tvd(samples, outcomes, probs)
    println("T3 TVD (Householder G=I) = ", round(tvd, digits=4))
    @test tvd < 0.05
end

# ---------------------------------------------------------------------------
# T4 — n=1 sanity. A single photon can't interfere with itself, so the
# output distribution is `|U_bsj[1, j]|²` in mode j (or `|U_phys[j, 1]|²`).
# The Gram matrix is a scalar 1×1 = [1], so Householder has NO internal
# DOF; this isolates the "build_block_diagonal_interferometers + shuffle"
# branch without any photon–photon interference.
# ---------------------------------------------------------------------------
@testset "T4 [diagnostic]: single-photon Householder matches |U|² row" begin
    Random.seed!(SEED)
    n, m, N = 1, 5, 5_000    # only 5 outcomes — small N is fine
    interf = RandHaar(m)
    input  = Input{UserDefinedGramMatrix}(first_modes(n, m),
                                          ComplexF64[1.0;;])

    outcomes, probs = exact_distribution(input, interf)
    samples = draw_householder(input, interf, N)
    tvd = empirical_tvd(samples, outcomes, probs)
    println("T4 TVD (Householder n=1) = ", round(tvd, digits=4))
    @test tvd < 0.05
end

# ---------------------------------------------------------------------------
# T5 — HOM dip. n=2, m=2, 50:50 BS, `S = [1 s; conj(s) 1]`:
#     P((1,1)) = (1 - |s|²) / 2
# Try s ∈ {0 (distinguishable), 0.5, 1 (bosonic)}. At s=0 we expect
# P(1,1) = 1/2; at s=1, P(1,1) = 0. Any convention flip shows up here
# because the beamsplitter is symmetric so U = transpose(U), meaning
# physics vs. BSJ convention cannot be the cause.
# ---------------------------------------------------------------------------
@testset "T5 [diagnostic]: HOM dip, n=2, m=2, varying overlap" begin
    Random.seed!(SEED)
    BS = ComplexF64[1.0  1.0;
                    1.0 -1.0] ./ sqrt(2)
    interf = UserDefinedInterferometer(BS)     # BS is symmetric, no transpose needed

    for s in (0.0, 0.5, 1.0)
        S = ComplexF64[1.0       s;
                       conj(s)   1.0]
        input = Input{UserDefinedGramMatrix}(first_modes(2, 2), S)
        outcomes, probs = exact_distribution(input, interf)

        N = 5_000    # only 3 outcomes (2,0),(1,1),(0,2) — small N is fine
        samples = draw_householder(input, interf, N)
        tvd = empirical_tvd(samples, outcomes, probs)
        p11_expected = (1 - abs2(s)) / 2
        idx_11 = findfirst(==([1, 1]), outcomes)
        @test probs[idx_11] ≈ p11_expected atol=1e-8   # exact formula sanity
        p11_empirical = count(==([1, 1]), samples) / N
        println("T5 s=$s  P(1,1) expected=$(round(p11_expected, digits=3))  ",
                "empirical=$(round(p11_empirical, digits=3))  TVD=$(round(tvd, digits=3))")
        @test tvd < 0.05
    end
end

# ---------------------------------------------------------------------------
# T6 — structural check on `build_householder_interferometer`. Whatever
# convention is intended, the output MUST be unitary; else the
# downstream Clifford call is meaningless. This is a pure linear-algebra
# check, no sampling involved.
# ---------------------------------------------------------------------------
@testset "T6 [structural]: expanded interferometer is unitary" begin
    Random.seed!(SEED)
    n, m = 3, 5
    U = RandHaar(m).U                        # m×m unitary
    σ = collect(1:n)                         # first_modes(n, m)
    # Try Gram matrices of different rank.
    for S in (ComplexF64.(ones(n, n)),
              ComplexF64.(Matrix(I, n, n)),
              rand_gram_matrix_from_orthonormal_basis(n, 2))
        C = BosonSampling.gram_to_coefficients(S)
        r = size(C, 2)
        full = BosonSampling.build_householder_interferometer(C, U, σ, n, m, r)
        @test size(full) == (r * m, r * m)   # DOF-major: r blocks of m modes
        dev = norm(full * full' - I)
        println("T6 rank r=$r  ‖V V† − I‖ = ", round(dev, sigdigits=3))
        @test dev < 1e-8
    end
end

# ---------------------------------------------------------------------------
# T7 — inspect the distribution on the EXPANDED (r·m) space. Householder
# runs Clifford on an `r·m`-mode Bosonic sampler (DOF-major layout:
# position (k, j) = (k-1)·m + j) and bins across DOF blocks. T7a checks
# that expanded Clifford matches expanded full_distribution; T7b checks
# that binning the exact expanded distribution reproduces the m-mode
# Bosonic distribution in the S = ones limit (r = 1 ⇒ expansion is a
# no-op, so this holds exactly).
# ---------------------------------------------------------------------------
@testset "T7 [diagnostic]: expanded-space distribution before binning" begin
    Random.seed!(SEED)
    n, m, N = 3, 5, 8_000
    U = RandHaar(m).U
    S = ComplexF64.(ones(n, n))              # Bosonic limit → r = 1
    C = BosonSampling.gram_to_coefficients(S)
    r = size(C, 2)
    σ = collect(1:n)                         # first_modes
    full = BosonSampling.build_householder_interferometer(C, U, σ, n, m, r)
    interf_expanded = UserDefinedInterferometer(full)

    # expanded input state (DOF-major): photon i at (k=1, j=σ[i]) = σ[i]
    occ = zeros(Int, r * m)
    for photon in 1:n
        occ[σ[photon]] = 1
    end
    input_expanded = Input{Bosonic}(ModeOccupation(occ))

    outcomes_exp, probs_exp = exact_distribution(input_expanded, interf_expanded)
    samples_exp = [cliffords_sampler(input=input_expanded, interf=interf_expanded)
                   for _ in 1:N]
    tvd_exp = empirical_tvd(samples_exp, outcomes_exp, probs_exp)
    println("T7 TVD (expanded Clifford vs expanded full_dist) = ", round(tvd_exp, digits=4))
    @test tvd_exp < 0.05

    # Bin each expanded outcome to m physical modes (sum across DOF blocks)
    # and compare to the exact m-mode Bosonic distribution.
    bin_to_m(v) = [sum(v[(k - 1) * m + j] for k in 1:r) for j in 1:m]
    outcomes_binned_counts = Dict{Vector{Int}, Float64}()
    for (k, o) in enumerate(outcomes_exp)
        b = bin_to_m(o)
        outcomes_binned_counts[b] = get(outcomes_binned_counts, b, 0.0) + probs_exp[k]
    end

    input_B = Input{Bosonic}(first_modes(n, m))
    out_B, probs_B = exact_distribution(input_B, UserDefinedInterferometer(U))

    tvd_binned = 0.0
    for (k, o) in enumerate(out_B)
        p_binned = get(outcomes_binned_counts, o, 0.0)
        tvd_binned += abs(p_binned - probs_B[k])
    end
    tvd_binned /= 2
    println("T7 TVD (expanded full_dist, binned, vs m-mode Bosonic) = ",
            round(tvd_binned, digits=4))
    @test tvd_binned < 1e-6              # pure linear-algebra equality — no stats
end

# ---------------------------------------------------------------------------
# T8 — `sample_multiple` matches `householder_sampler` called N times.
# If T8 fails, the bug is in the batch wrapper; if T8 passes but T2
# fails, the bug is in `householder_sampler` itself.
# ---------------------------------------------------------------------------
@testset "T8 [diagnostic]: sample_multiple agrees with looped householder" begin
    Random.seed!(SEED)
    n, m, N = 3, 5, 8_000
    interf = RandHaar(m)
    input  = Input{UserDefinedGramMatrix}(first_modes(n, m),
                                          ComplexF64.(ones(n, n)))

    samples_loop  = draw_householder(input, interf, N)
    samples_batch = sample_multiple(input, interf, N)
    # same distribution, not same samples (they use different RNG draws)
    outcomes_loop  = sort(collect(keys(histogram(samples_loop))))
    outcomes_batch = sort(collect(keys(histogram(samples_batch))))
    @test outcomes_loop == outcomes_batch

    h_loop  = histogram(samples_loop)
    h_batch = histogram(samples_batch)
    all_outcomes = union(keys(h_loop), keys(h_batch))
    tvd = 0.0
    for o in all_outcomes
        tvd += abs(get(h_loop, o, 0) / N - get(h_batch, o, 0) / N)
    end
    tvd /= 2
    println("T8 TVD (householder_sampler loop vs sample_multiple) = ",
            round(tvd, digits=4))
    @test tvd < 0.05
end

end  # outer @testset "Householder sampler: distribution correctness"
