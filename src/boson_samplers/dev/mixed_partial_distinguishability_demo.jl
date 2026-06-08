"""
Validation: Mixed-State Partial Distinguishability Sampler

Two checks for `mixed_partial_distinguishability_sampler` /
`sample_mixed_multiple`:

  A. Rank-1 reduction — when each ρ_i = |ψ_i⟩⟨ψ_i| is pure, the mixed sampler
     must reproduce the pure Householder sampler (same Gram matrix). Compared
     by total variation distance (TVD) of the two empirical histograms.

  B. Exact mixture — for a small system with genuinely mixed ρ_i, the empirical
     histogram must match the exact mixed distribution computed as the convex
     combination Σ (∏λ) · p_pure over all eigenstate combinations, with each
     p_pure obtained from `compute_probability!`.

Run with:  julia --project -e 'include("src/boson_samplers/dev/mixed_partial_distinguishability_demo.jl")'
"""

using BosonSampling
using Random
using StatsBase
using LinearAlgebra

Random.seed!(1234)

# ---- helpers ---------------------------------------------------------------

"Empirical probability vector over `outputs` (vector of ModeOccupation) from a
list of sampled state vectors."
function empirical_distribution(samples::Vector{Vector{Int}}, outputs)
    counts = zeros(Int, length(outputs))
    index = Dict(out.state => k for (k, out) in enumerate(outputs))
    for s in samples
        haskey(index, s) && (counts[index[s]] += 1)
    end
    return counts ./ sum(counts)
end

tvd(p, q) = 0.5 * sum(abs.(p .- q))

"Random normalised complex vector in C^d."
rand_pure(d) = (v = randn(ComplexF64, d); v / norm(v))

"Exact mixed-state distribution as a convex combination of pure partially
distinguishable probabilities, enumerated over all eigenstate choices."
function exact_mixed_distribution(prep, r::ModeOccupation, interf, outputs)
    n = length(prep)
    probs = zeros(Float64, length(outputs))
    ranges = [1:length(prep[i].probs) for i in 1:n]
    for combo in Iterators.product(ranges...)
        weight = prod(prep[i].probs[combo[i]] for i in 1:n)
        weight == 0 && continue
        W = reduce(hcat, [prep[i].vecs[:, combo[i]] for i in 1:n])
        S = gram_from_vectors(W)
        input = Input{UserDefinedGramMatrix}(r, S)
        for (k, out) in enumerate(outputs)
            ev = Event(input, FockDetection(out), interf)
            compute_probability!(ev)
            probs[k] += weight * real(ev.proba_params.probability)
        end
    end
    return probs
end

# ---- Test A: rank-1 reduces to the pure Householder sampler -----------------

println("=" ^ 60)
println("TEST A — rank-1 ρ_i must match the pure Householder sampler")
println("=" ^ 60)

let
    n, m, d = 3, 5, 3
    r = first_modes(n, m)
    interf = RandHaar(m)
    N = 200_000

    ψ = [rand_pure(d) for _ in 1:n]
    ρ_list = [ψ[i] * ψ[i]' for i in 1:n]          # rank-1 density matrices
    W = reduce(hcat, ψ)
    S_pure = gram_from_vectors(W)

    outputs = ModeOccupation.(all_mode_configurations(n, m; only_photon_number_conserving = true))

    samples_mixed = sample_mixed_multiple(ρ_list, r, interf, N)
    samples_pure  = sample_householder_multiple(Input{UserDefinedGramMatrix}(r, S_pure), interf, N)

    p_mixed = empirical_distribution(samples_mixed, outputs)
    p_pure  = empirical_distribution(samples_pure, outputs)

    d_tv = tvd(p_mixed, p_pure)
    println("n=$n, m=$m, d=$d, samples=$N")
    println("TVD(mixed-sampler, pure-sampler) = $(round(d_tv, digits=4))")
    println(d_tv < 0.02 ? "✅ PASS (rank-1 reduces to pure)" : "❌ FAIL")
end

# ---- Test B: genuine mixture matches the exact convex combination -----------

println()
println("=" ^ 60)
println("TEST B — mixed ρ_i must match the exact convex-combination law")
println("=" ^ 60)

let
    n, m, d = 2, 3, 2
    r = first_modes(n, m)
    interf = RandHaar(m)
    N = 300_000

    # Per-photon mixed states in a shared 2-dim internal basis.
    # Photon 1: biased mixture of |0> and |1>; photon 2: mixture of |+> and |->.
    e0 = ComplexF64[1, 0]; e1 = ComplexF64[0, 1]
    plus = (e0 + e1) / sqrt(2); minus = (e0 - e1) / sqrt(2)
    ρ1 = 0.7 * (e0 * e0') + 0.3 * (e1 * e1')
    ρ2 = 0.6 * (plus * plus') + 0.4 * (minus * minus')
    ρ_list = [ρ1, ρ2]

    prep = eigendecompose_density_matrices(ρ_list)
    outputs = ModeOccupation.(all_mode_configurations(n, m; only_photon_number_conserving = true))

    p_exact = exact_mixed_distribution(prep, r, interf, outputs)
    samples = sample_mixed_multiple(ρ_list, r, interf, N)
    p_emp = empirical_distribution(samples, outputs)

    d_tv = tvd(p_emp, p_exact)
    println("n=$n, m=$m, d=$d, samples=$N")
    println("Σ p_exact = $(round(sum(p_exact), digits=6)) (should be 1)")
    for (k, out) in enumerate(outputs)
        println("  $(out.state): empirical $(round(p_emp[k], digits=4))  exact $(round(p_exact[k], digits=4))")
    end
    println("TVD(empirical, exact) = $(round(d_tv, digits=4))")
    println(d_tv < 0.01 ? "✅ PASS (matches exact mixed law)" : "❌ FAIL")
end
