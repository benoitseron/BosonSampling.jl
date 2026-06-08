"""
Mixed-State Partial Distinguishability Sampler

Extends the pure-state partial distinguishability sampler
([`householder_sampler`](@ref)) to photons whose internal (distinguishability)
degrees of freedom are described by *density matrices* ρ_i rather than pure
states |ψ_i⟩.

Key observation: the boson sampling output distribution is linear in each ρ_i,
so a mixed internal state is operationally a classical mixture over pure
configurations. Writing each density matrix in its eigenbasis

    ρ_i = Σ_α λ_α^(i) |v_α^(i)⟩⟨v_α^(i)|,   λ_α^(i) ≥ 0, Σ_α λ_α^(i) = 1,

the joint state ρ = ⊗_i ρ_i expands as

    ρ = Σ_{α_1,...,α_n} (∏_i λ_{α_i}^(i)) ⊗_i |v_{α_i}^(i)⟩⟨v_{α_i}^(i)|,

and consequently

    p(s | U, ρ) = Σ_{α_1,...,α_n} (∏_i λ_{α_i}^(i)) · p(s | U, {|v_{α_i}^(i)⟩}),

where each term on the right is an ordinary *pure* partially distinguishable
boson sampling probability already handled by [`householder_sampler`](@ref).

Exact sampling algorithm (one sample):

1. (Preprocessing, once) Diagonalise each ρ_i → eigenvalues {λ_α^(i)} and
   eigenvectors {|v_α^(i)⟩}. Cost O(n d³), d = internal dimension.
2. For each photon i independently, draw α_i ∼ Categorical(λ_·^(i)). This
   assigns photon i the pure internal state |v_{α_i}^(i)⟩.
3. Build the Gram matrix S of the drawn pure states and feed it to the
   existing pure-state Householder sampler.

The output of each call is one *exact* sample from the mixed-state
distribution; the only overhead beyond the pure sampler is the O(n)
categorical draws and an O(n²d) Gram matrix build per sample.

All density matrices must live in a shared internal Hilbert space of the same
dimension d so that the overlaps ⟨v|w⟩ are well defined. Rank-1 ρ_i (pure
states) reduce exactly to the pure-state path.
"""

using LinearAlgebra
using StatsBase

"""
    gram_from_vectors(W::AbstractMatrix) -> Matrix

Gram matrix of the columns of `W`. With `W` a `d×n` matrix whose column `i`
is the internal-state vector |v_i⟩, returns the `n×n` Hermitian matrix
`S[i,j] = ⟨v_i|v_j⟩ = Σ_k conj(W[k,i]) W[k,j]`, i.e. `S = W' * W`. This is the
physics-convention Gram matrix consumed by [`householder_sampler`](@ref).
"""
gram_from_vectors(W::AbstractMatrix) = W' * W

"""
    eigendecompose_density_matrices(ρ_list; atol=ATOL)
        -> Vector{NamedTuple{(:probs, :vecs)}}

Diagonalise each per-photon density matrix in `ρ_list` (a length-`n` vector of
`d×d` Hermitian, PSD, trace-1 matrices in a shared internal basis). Returns,
for each photon, a `NamedTuple` with

  - `probs::Vector{Float64}` — eigenvalues (mixture weights), clamped to ≥ 0
    and renormalised to sum to 1,
  - `vecs::Matrix{ComplexF64}` — corresponding eigenvectors as columns.

Validates that every matrix is Hermitian, square, of common dimension `d`,
has non-negative eigenvalues (up to `atol`) and unit trace (up to `atol`).
"""
function eigendecompose_density_matrices(ρ_list::AbstractVector{<:AbstractMatrix}; atol::Real = ATOL)

    n = length(ρ_list)
    n == 0 && error("ρ_list must contain at least one density matrix.")

    d = size(ρ_list[1], 1)

    prep = Vector{NamedTuple{(:probs, :vecs), Tuple{Vector{Float64}, Matrix{ComplexF64}}}}(undef, n)

    for (i, ρ) in enumerate(ρ_list)
        size(ρ, 1) == size(ρ, 2) || error("ρ[$i] must be square.")
        size(ρ, 1) == d || error("all density matrices must share the same internal dimension d=$d (ρ[$i] has dimension $(size(ρ,1))).")
        ishermitian(ρ) || isapprox(ρ, ρ', atol = atol) || error("ρ[$i] must be Hermitian.")
        isapprox(real(tr(ρ)), 1.0, atol = atol) || error("ρ[$i] must have unit trace (got $(tr(ρ))).")

        F = eigen(Hermitian(Matrix{ComplexF64}(ρ)))
        λ = real.(F.values)
        all(λ .>= -atol) || error("ρ[$i] is not positive semi-definite (min eigenvalue $(minimum(λ))).")

        λ = max.(λ, 0.0)
        s = sum(λ)
        s > 0 || error("ρ[$i] has zero trace after clamping.")
        λ ./= s

        prep[i] = (probs = λ, vecs = Matrix{ComplexF64}(F.vectors))
    end

    return prep
end

"""
    draw_gram_matrix(prep) -> Matrix

Draw one mixture component: for each photon, sample an eigenstate index
`α_i ∼ Categorical(prep[i].probs)`, collect the chosen eigenvectors as the
columns of a `d×n` matrix `W`, and return the Gram matrix `S = W' * W` of the
drawn pure internal states.
"""
function draw_gram_matrix(prep::AbstractVector)
    n = length(prep)
    d = size(prep[1].vecs, 1)
    W = Matrix{ComplexF64}(undef, d, n)
    for i in 1:n
        α = StatsBase.sample(Weights(prep[i].probs))
        W[:, i] = prep[i].vecs[:, α]
    end
    return gram_from_vectors(W)
end

"""
    mixed_partial_distinguishability_sampler(ρ_list, r::ModeOccupation,
                                             interf::Interferometer; prep=nothing)
        -> ModeOccupation

Draw one exact sample from a boson sampling experiment with `n` photons whose
internal degrees of freedom are described by the per-photon density matrices
`ρ_list` (length `n`, each `d×d`), entering the interferometer `interf` at the
input occupation `r` (distinct occupied modes, no input bunching).

If `prep` (the output of [`eigendecompose_density_matrices`](@ref)) is supplied
the eigendecomposition is reused; otherwise it is computed on the fly. For many
samples with the same `ρ_list`, prefer [`sample_mixed_multiple`](@ref).
"""
function mixed_partial_distinguishability_sampler(ρ_list::AbstractVector{<:AbstractMatrix},
                                                  r::ModeOccupation,
                                                  interf::Interferometer;
                                                  prep = nothing)

    length(ρ_list) == r.n || error("number of density matrices ($(length(ρ_list))) must equal photon number ($(r.n)).")
    prep === nothing && (prep = eigendecompose_density_matrices(ρ_list))

    S = draw_gram_matrix(prep)
    input = Input{UserDefinedGramMatrix}(r, S)
    ev = Event(input, FockSample(), interf)

    return householder_sampler(ev)
end

"""
    mixed_partial_distinguishability_sampler(ev::Event{MixedDensityMatrices, FockSample})
        -> ModeOccupation

Event-based entry point used by `sample!`. Reads the per-photon density
matrices from the input's Gram matrix and draws one exact sample.
"""
function mixed_partial_distinguishability_sampler(ev::Event{MixedDensityMatrices, TOut}) where {TOut <: FockSample}
    input = ev.input_state
    return mixed_partial_distinguishability_sampler(input.G.density_matrices, input.r, ev.interferometer)
end

"""
    sample_mixed_multiple(ρ_list, r::ModeOccupation, interf::Interferometer,
                          n_samples::Int; threaded=false, show_progress=false)
        -> Vector{Vector{Int}}

Generate `n_samples` exact samples for a fixed mixed-state setup. The
eigendecomposition of the `ρ_i` is computed once and reused; the Gram matrix
(and hence the expanded interferometer) is redrawn on every sample, since each
mixture component generally has a different Gram matrix.
"""
function sample_mixed_multiple(ρ_list::AbstractVector{<:AbstractMatrix},
                               r::ModeOccupation,
                               interf::Interferometer,
                               n_samples::Int;
                               threaded::Bool = false,
                               show_progress::Bool = false)

    prep = eigendecompose_density_matrices(ρ_list)
    samples = Vector{Vector{Int}}(undef, n_samples)

    if threaded
        Threads.@threads for i in 1:n_samples
            samples[i] = mixed_partial_distinguishability_sampler(ρ_list, r, interf; prep = prep).state
        end
    else
        iter = show_progress ? ProgressBar(1:n_samples) : (1:n_samples)
        for i in iter
            samples[i] = mixed_partial_distinguishability_sampler(ρ_list, r, interf; prep = prep).state
        end
    end

    return samples
end
