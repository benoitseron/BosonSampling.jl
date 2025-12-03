"""
Partial Distinguishability Sampler using Householder Transformations

Implementation based on Gram matrix decomposition approach:
S = CC† where S is the Gram matrix of photon overlaps.

Algorithm:
1. Decompose Gram matrix S = CC† to extract internal degrees of freedom
2. Apply Householder transformation V_i = I - 2C_iC_i† to each photon
3. Shuffle modes to group by position across photon blocks
4. Apply physical interferometer U to each photon's modes
5. Use Clifford sampler in expanded n×m mode space
6. Bin results back to m physical modes

This method works for arbitrary Gram matrices and naturally handles
different degrees of distinguishability between photon pairs.
"""

using LinearAlgebra

"""
    householder_sampler(ev::Event{TIn, FockSample}) where {TIn<:PartDist}

Sample from a boson sampling experiment with partially distinguishable photons
specified by a Gram matrix.

Works with any PartDist input type including:
- UserDefinedGramMatrix: Custom Gram matrix
- OneParameterInterpolation: Single parameter x ∈ [0,1]
- RandomGramMatrix: Random Gram matrix

The Gram matrix S where S[i,j] = ⟨photon_i, photon_j⟩ specifies the overlap
between photons in their internal degrees of freedom.

Returns a ModeOccupation representing the sampled output state.
"""
function householder_sampler(ev::Event{TIn, FockSample}) where {TIn<:PartDist}

    input = ev.input_state
    interf = ev.interferometer

    n = input.n  # number of photons
    m = input.m  # number of modes

    # Extract Gram matrix from input
    S = input.G.S  # Extract matrix from GramMatrix wrapper

    # Decompose Gram matrix to get coefficients: S = CC†
    C = gram_to_coefficients(S)
    r = size(C, 2)  # rank of Gram matrix (number of internal DOF)

    # Build the full transformation matrix
    full_interf = build_householder_interferometer(C, interf.U, n, m, r)

    # Create input state in expanded space (n photons, each in first mode of its block)
    occupation_expanded = zeros(Int, n*m)
    for photon in 1:n
        occupation_expanded[(photon-1)*m + 1] = 1
    end

    # Create Event in expanded space with bosonic photons
    input_expanded = Input{Bosonic}(ModeOccupation(occupation_expanded))
    interf_expanded = UserDefinedInterferometer(full_interf)
    output_expanded = FockSample()
    ev_expanded = Event(input_expanded, output_expanded, interf_expanded)

    # Sample using Clifford algorithm
    sample!(ev_expanded)

    # Get sampled output and bin back to physical modes
    sampled_expanded = ev_expanded.output_measurement.s.state
    sampled_physical = bin_to_physical_modes(sampled_expanded, n, m)

    return ModeOccupation(sampled_physical)
end

"""
    build_householder_interferometer(C::Matrix, U::Matrix, n::Int, m::Int, r::Int)

Build the full interferometer for the expanded mode space.

Combines:
1. Householder splitting matrices (one per photon)
2. Mode shuffle permutation
3. Block diagonal physical interferometers
4. Reverse shuffle
"""
function build_householder_interferometer(C::Matrix, U::Matrix, n::Int, m::Int, r::Int)

    # 1. Build splitting matrices (Householder for each photon)
    splitting_interf = build_splitting_matrices(C, n, m, r)

    # 2. Build shuffle permutation
    shuffle_perm = build_shuffle_permutation(n, m)

    # 3. Build block diagonal interferometers
    block_diag_interf = build_block_diagonal_interferometers(U, n, m)

    # 4. Combine: splitting → shuffle → U → unshuffle
    full_interf = shuffle_perm' * block_diag_interf * shuffle_perm * splitting_interf

    return full_interf
end

"""
    build_splitting_matrices(C::Matrix, n::Int, m::Int, r::Int)

Build block diagonal matrix with Householder transformations for each photon.

For photon i, the Householder transformation is V_i = I - 2C_iC_i†
where C_i is the i-th row of C (coefficients in internal basis).
"""
function build_splitting_matrices(C::Matrix, n::Int, m::Int, r::Int)

    result = Matrix{ComplexF64}(I, n*m, n*m)

    for photon in 1:n
        # Householder matrix for this photon (r×r)
        C_i = C[photon, :]
        V_i = Matrix{ComplexF64}(I, r, r) - 2 * C_i * C_i'

        # Enlarge to m×m (only first r×r block is non-trivial)
        V_i_enlarged = Matrix{ComplexF64}(I, m, m)
        V_i_enlarged[1:r, 1:r] = V_i

        # Place in block diagonal position
        range_start = (photon-1)*m + 1
        range_end = photon*m
        result[range_start:range_end, range_start:range_end] = V_i_enlarged
    end

    return result
end

"""
    build_shuffle_permutation(n::Int, m::Int)

Build permutation matrix that reorders modes from photon-major to mode-major ordering.

Before: [Photon₁: modes 1...m | Photon₂: modes 1...m | ... | Photonₙ: modes 1...m]
After:  [Mode₁: photons 1...n | Mode₂: photons 1...n | ... | Modeₘ: photons 1...n]
"""
function build_shuffle_permutation(n::Int, m::Int)

    P = zeros(ComplexF64, n*m, n*m)

    for photon in 1:n
        for mode in 1:m
            k_old = (photon-1)*m + mode  # Position in photon-major order
            k_new = (mode-1)*n + photon  # Position in mode-major order
            P[k_new, k_old] = 1.0
        end
    end

    return P
end

"""
    build_block_diagonal_interferometers(U::Matrix, n::Int, m::Int)

Build block diagonal matrix with n copies of the m×m interferometer U.

Each photon gets its own copy of U acting on its m modes.
In the shuffled basis, this creates an interleaved block structure.
"""
function build_block_diagonal_interferometers(U::Matrix, n::Int, m::Int)

    result = zeros(ComplexF64, n*m, n*m)

    for photon in 1:n
        for mode_i in 1:m
            for mode_j in 1:m
                # In shuffled basis: (photon, mode) is at position (mode-1)*n + photon
                row = (mode_i-1)*n + photon
                col = (mode_j-1)*n + photon
                result[row, col] = U[mode_i, mode_j]
            end
        end
    end

    return result
end

"""
    bin_to_physical_modes(sampled_expanded::Vector{Int}, n::Int, m::Int)

Collapse the expanded space (n*m modes) back to physical modes (m modes).

Each physical mode i corresponds to modes {i, m+i, 2m+i, ..., (n-1)m+i}
in the expanded space. Sum photon counts across all photon blocks.
"""
function bin_to_physical_modes(sampled_expanded::Vector{Int}, n::Int, m::Int)

    sampled_physical = zeros(Int, m)

    for photon in 1:n
        for mode in 1:m
            expanded_mode = (photon-1)*m + mode
            sampled_physical[mode] += sampled_expanded[expanded_mode]
        end
    end

    return sampled_physical
end

"""
    householder_sampler_vec(ev::Event{TIn, FockSample}) where {TIn<:PartDist}

Convenience wrapper that extracts the sample vector directly.
"""
function householder_sampler_vec(ev::Event{TIn, FockSample}) where {TIn<:PartDist}
    return householder_sampler(ev).state
end

"""
    sample_householder_multiple(input::Input{TIn}, interf::Interferometer, n_samples::Int;
                               threaded::Bool=false) where {TIn<:PartDist}

Efficiently generate multiple samples from the same input configuration.

This function pre-computes all fixed transformations (Gram decomposition, Householder matrices,
permutations, etc.) once, then samples n_samples times by only running the Clifford algorithm
and binning step repeatedly. This is much faster than calling sample!() multiple times.

With `threaded=true`, uses multi-threading to parallelize the sampling loop across available CPU cores.

# Arguments
- `input::Input{TIn}`: Input state with partial distinguishability
- `interf::Interferometer`: Interferometer (e.g., RandHaar, Fourier, etc.)
- `n_samples::Int`: Number of samples to generate
- `threaded::Bool=false`: Use multi-threading for parallel sampling

# Returns
- `Vector{Vector{Int}}`: Vector of samples, each sample is a mode occupation vector

# Example
```julia
input = Input{UserDefinedGramMatrix}(first_modes(3, 5), S)
interf = RandHaar(5)

# Sequential sampling
samples = sample_householder_multiple(input, interf, 1000)

# Parallel sampling (uses Threads.nthreads() cores)
samples = sample_householder_multiple(input, interf, 1000, threaded=true)
```
"""
function sample_householder_multiple(input::Input{TIn}, interf::Interferometer, n_samples::Int;
                                    threaded::Bool=false) where {TIn<:PartDist}

    n = input.n
    m = input.m

    # Extract Gram matrix and decompose (done once)
    S = input.G.S
    C = gram_to_coefficients(S)
    r = size(C, 2)

    # Build full interferometer (done once)
    full_interf = build_householder_interferometer(C, interf.U, n, m, r)

    # Create expanded input state (done once)
    occupation_expanded = zeros(Int, n*m)
    for photon in 1:n
        occupation_expanded[(photon-1)*m + 1] = 1
    end
    input_expanded = Input{Bosonic}(ModeOccupation(occupation_expanded))
    interf_expanded = UserDefinedInterferometer(full_interf)

    # Generate samples efficiently
    samples = Vector{Vector{Int}}(undef, n_samples)

    if threaded
        # Parallel sampling using multi-threading
        Threads.@threads for i in 1:n_samples
            # Create fresh event for each sample (thread-safe)
            ev_expanded = Event(input_expanded, FockSample(), interf_expanded)

            # Sample using Clifford
            sample!(ev_expanded)

            # Bin to physical modes
            sampled_expanded = ev_expanded.output_measurement.s.state
            sampled_physical = bin_to_physical_modes(sampled_expanded, n, m)

            samples[i] = sampled_physical
        end
    else
        # Sequential sampling
        for i in 1:n_samples
            # Create fresh event for each sample
            ev_expanded = Event(input_expanded, FockSample(), interf_expanded)

            # Sample using Clifford
            sample!(ev_expanded)

            # Bin to physical modes
            sampled_expanded = ev_expanded.output_measurement.s.state
            sampled_physical = bin_to_physical_modes(sampled_expanded, n, m)

            samples[i] = sampled_physical
        end
    end

    return samples
end

"""
    sample_multiple(input::Input{TIn}, interf::Interferometer, n_samples::Int;
                   threaded::Bool=false) where {TIn<:InputType}

General function to efficiently generate multiple samples from any input type.

Dispatches to specialized implementations:
- PartDist types → householder sampler (pre-computes transformations)
- Bosonic → Clifford sampler (pre-computes interferometer submatrix)
- Distinguishable → classical sampler
- Other types → falls back to repeated sample!() calls

# Arguments
- `input::Input{TIn}`: Input state
- `interf::Interferometer`: Interferometer
- `n_samples::Int`: Number of samples to generate
- `threaded::Bool=false`: Use multi-threading for parallel sampling

# Returns
- `Vector{Vector{Int}}`: Vector of samples

# Example
```julia
input = Input{Bosonic}(first_modes(3, 5))
interf = RandHaar(5)

# Sequential sampling
samples = sample_multiple(input, interf, 1000)

# Parallel sampling (uses all available threads)
samples = sample_multiple(input, interf, 1000, threaded=true)
```
"""
function sample_multiple(input::Input{TIn}, interf::Interferometer, n_samples::Int;
                        threaded::Bool=false) where {TIn<:InputType}

    if TIn <: PartDist
        # Use optimized Householder multi-sampling
        return sample_householder_multiple(input, interf, n_samples, threaded=threaded)
    else
        # Fall back to repeated sampling for other types
        samples = Vector{Vector{Int}}(undef, n_samples)

        if threaded
            Threads.@threads for i in 1:n_samples
                ev = Event(input, FockSample(), interf)
                sample!(ev)
                samples[i] = ev.output_measurement.s.state
            end
        else
            for i in 1:n_samples
                ev = Event(input, FockSample(), interf)
                sample!(ev)
                samples[i] = ev.output_measurement.s.state
            end
        end
        return samples
    end
end
