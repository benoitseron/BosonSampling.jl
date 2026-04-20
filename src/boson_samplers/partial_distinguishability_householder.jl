"""
Partial Distinguishability Sampler using Householder Transformations

Implementation based on the Gram matrix decomposition approach:
S = C C* where S[i, j] is the overlap between photons i and j in their
internal degrees of freedom. C is n×r, one row per photon, r = rank(S).

Algorithm (reducing to Clifford on an expanded Bosonic state):

1. Decompose Gram matrix S = C C* (`gram_to_coefficients`).
2. Expand the mode space to r·m modes, organised as r DOF blocks of m
   spatial modes each. Layout is DOF-major: position (k, j) = (k-1)·m + j,
   with k ∈ 1..r (DOF) and j ∈ 1..m (spatial mode).
3. Place photon i at (k=1, j=σ(i)), i.e. the DOF=1 copy of its physical
   input spatial mode σ(i).
4. Splitting: for each photon i, apply a per-photon r×r unitary V_i on
   the DOF coordinate at spatial mode σ(i). V_i is chosen so that its
   first column equals C[i, :] (i.e. V_i @ e_1 = C[i, :] in physics
   convention). Because the splitting matrix sits inside a package-
   convention interferometer (which Clifford transposes), the sub-block
   stored is `transpose(V_i)`.
5. Apply r block-diagonal copies of the physical interferometer U, one
   per DOF block. Photons within the same DOF block interfere bosonically;
   photons in disjoint DOF blocks do not.
6. Sample via Clifford on the expanded Bosonic state, then bin the
   expanded outcome back to m physical modes by summing across DOF.

Reduces correctly to Clifford Bosonic when r = 1 (e.g. S = ones) and to
the classical Distinguishable sampler when r = n (e.g. S = I).
"""

using LinearAlgebra

"""
    householder_sampler(ev::Event{TIn, FockSample}) where {TIn<:PartDist}

Sample from a boson sampling experiment with partially distinguishable
photons specified by a Gram matrix.

Requires the input to have distinct occupied spatial modes (no bunching
at input). Returns a `ModeOccupation` of the sampled output state on m
physical modes.
"""
function householder_sampler(ev::Event{TIn, FockSample}) where {TIn<:PartDist}

    input = ev.input_state
    interf = ev.interferometer

    n = input.n
    m = input.m

    # Photon i is at spatial mode σ[i] (1-based); σ is the sequence of
    # occupied input modes (same indexing as the Gram matrix rows).
    σ = fill_arrangement(input)
    @argcheck length(σ) == n
    @argcheck length(unique(σ)) == n "Householder sampler requires distinct input modes (no bunched inputs)."

    # Gram matrix S = C C* → C is n × r, r = rank(S)
    S = input.G.S
    C = gram_to_coefficients(S)
    r = size(C, 2)

    full_interf = build_householder_interferometer(C, interf.U, σ, n, m, r)

    occupation_expanded = zeros(Int, r * m)
    for photon in 1:n
        occupation_expanded[σ[photon]] = 1        # position (k=1, j=σ[photon])
    end

    input_expanded = Input{Bosonic}(ModeOccupation(occupation_expanded))
    interf_expanded = UserDefinedInterferometer(full_interf)
    ev_expanded = Event(input_expanded, FockSample(), interf_expanded)

    sample!(ev_expanded)

    sampled_expanded = ev_expanded.output_measurement.s.state
    sampled_physical = bin_to_physical_modes(sampled_expanded, r, m)

    return ModeOccupation(sampled_physical)
end

"""
    build_householder_interferometer(C, U, σ, n, m, r)

Build the package-convention expanded interferometer (splitting · block_diag).
"""
function build_householder_interferometer(C::Matrix, U::Matrix,
                                          σ::AbstractVector{<:Integer},
                                          n::Int, m::Int, r::Int)
    splitting  = build_splitting_matrices(C, σ, n, m, r)
    block_diag = build_block_diagonal_interferometers(U, m, r)
    # Package convention: left-factor applied first. Splitting spreads each
    # photon across DOF blocks, then block_diag applies one U per block.
    return splitting * block_diag
end

"""
    build_splitting_matrices(C, σ, n, m, r)

Block-diagonal splitting unitary on the r·m expanded space. For each
photon i at spatial mode σ[i], acts as V_i on the DOF coordinate at
spatial mode σ[i] (positions {(k-1)·m + σ[i] : k=1..r}); identity on
unoccupied spatial modes.
"""
function build_splitting_matrices(C::Matrix, σ::AbstractVector{<:Integer},
                                  n::Int, m::Int, r::Int)
    result = Matrix{ComplexF64}(I, r * m, r * m)
    for photon in 1:n
        j = σ[photon]
        V_i = unitary_with_first_column(ComplexF64.(C[photon, :]))
        positions = [(k - 1) * m + j for k in 1:r]
        # Stored transposed: the surrounding matrix is package-convention
        # (M[in, out]), so the sub-block's first row must equal C[i, :].
        result[positions, positions] = transpose(V_i)
    end
    return result
end

"""
    build_block_diagonal_interferometers(U, m, r)

Build r copies of U on the diagonal of an (r·m)×(r·m) matrix (one per DOF
block). Columns 1..m = block 1, m+1..2m = block 2, etc.
"""
function build_block_diagonal_interferometers(U::Matrix, m::Int, r::Int)
    result = zeros(ComplexF64, r * m, r * m)
    for k in 1:r
        rng = (k - 1) * m + 1 : k * m
        result[rng, rng] = U
    end
    return result
end

"""
    unitary_with_first_column(c)

Return an r×r unitary V with V[:, 1] = c (expects ||c|| ≈ 1). Built via
QR of M = [c | e_{k≠i_drop}] and a phase correction so V[:, 1] equals c
exactly. Only the first column is load-bearing for the sampler; the
remaining columns just provide a unitary completion.
"""
function unitary_with_first_column(c::AbstractVector{ComplexF64})
    r = length(c)
    if r == 1
        return reshape(ComplexF64[c[1]], 1, 1)
    end
    i_drop = argmax(abs.(c))
    M = zeros(ComplexF64, r, r)
    M[:, 1] = c
    j = 1
    for k in 1:r
        k == i_drop && continue
        j += 1
        M[k, j] = 1.0 + 0.0im
    end
    Q, _ = qr(M)
    V = Matrix{ComplexF64}(Q)
    phase = c[i_drop] / V[i_drop, 1]
    V[:, 1] .*= phase
    return V
end

"""
    bin_to_physical_modes(sampled_expanded, r, m)

Collapse the r·m expanded occupation vector back to m physical modes by
summing across DOF blocks: physical mode j = Σ_k expanded[(k-1)·m + j].
"""
function bin_to_physical_modes(sampled_expanded::Vector{Int}, r::Int, m::Int)
    sampled_physical = zeros(Int, m)
    for k in 1:r, j in 1:m
        sampled_physical[j] += sampled_expanded[(k - 1) * m + j]
    end
    return sampled_physical
end

"""
    householder_sampler_vec(ev) -> Vector{Int}

Convenience wrapper that returns the sampled mode-occupation vector.
"""
function householder_sampler_vec(ev::Event{TIn, FockSample}) where {TIn<:PartDist}
    return householder_sampler(ev).state
end

"""
    sample_householder_multiple(input, interf, n_samples; threaded=false, show_progress=false)

Efficiently generate `n_samples` samples for a fixed `(input, interf)`.
Pre-computes the Gram decomposition, splitting, and block-diagonal U
once, then repeatedly calls Clifford on the expanded state.
"""
function sample_householder_multiple(input::Input{TIn}, interf::Interferometer, n_samples::Int;
                                     threaded::Bool=false, show_progress::Bool=false) where {TIn<:PartDist}
    n = input.n
    m = input.m

    σ = fill_arrangement(input)
    @argcheck length(σ) == n
    @argcheck length(unique(σ)) == n "Householder sampler requires distinct input modes (no bunched inputs)."

    S = input.G.S
    C = gram_to_coefficients(S)
    r = size(C, 2)

    full_interf = build_householder_interferometer(C, interf.U, σ, n, m, r)

    occupation_expanded = zeros(Int, r * m)
    for photon in 1:n
        occupation_expanded[σ[photon]] = 1
    end
    input_expanded = Input{Bosonic}(ModeOccupation(occupation_expanded))
    interf_expanded = UserDefinedInterferometer(full_interf)

    samples = Vector{Vector{Int}}(undef, n_samples)

    if threaded
        Threads.@threads for i in 1:n_samples
            ev_expanded = Event(input_expanded, FockSample(), interf_expanded)
            sample!(ev_expanded)
            samples[i] = bin_to_physical_modes(ev_expanded.output_measurement.s.state, r, m)
        end
    else
        iter = show_progress ? ProgressBar(1:n_samples) : (1:n_samples)
        for i in iter
            ev_expanded = Event(input_expanded, FockSample(), interf_expanded)
            sample!(ev_expanded)
            samples[i] = bin_to_physical_modes(ev_expanded.output_measurement.s.state, r, m)
        end
    end

    return samples
end

"""
    sample_multiple(input, interf, n_samples; threaded=false, show_progress=false)

Dispatches to specialised batch samplers:
- PartDist inputs → `sample_householder_multiple` (pre-computes expansion).
- Other input types → repeated `sample!()` calls.
"""
function sample_multiple(input::Input{TIn}, interf::Interferometer, n_samples::Int;
                         threaded::Bool=false, show_progress::Bool=false) where {TIn<:InputType}
    if TIn <: PartDist
        return sample_householder_multiple(input, interf, n_samples,
                                           threaded=threaded, show_progress=show_progress)
    end

    samples = Vector{Vector{Int}}(undef, n_samples)
    if threaded
        Threads.@threads for i in 1:n_samples
            ev = Event(input, FockSample(), interf)
            sample!(ev)
            samples[i] = ev.output_measurement.s.state
        end
    else
        iter = show_progress ? ProgressBar(1:n_samples) : (1:n_samples)
        for i in iter
            ev = Event(input, FockSample(), interf)
            sample!(ev)
            samples[i] = ev.output_measurement.s.state
        end
    end
    return samples
end
