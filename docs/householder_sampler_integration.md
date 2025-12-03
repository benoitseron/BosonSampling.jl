# Householder-Based Partial Distinguishability Sampler

## Overview

This implementation provides an integrated sampler for partially distinguishable photons using Householder transformations and Gram matrix decomposition. It works seamlessly with the existing `sample!()` interface for `Input{UserDefinedGramMatrix}` types.

## Algorithm

The algorithm is based on decomposing the Gram matrix S that specifies photon overlaps:

1. **Gram Matrix Decomposition**: S = CC† where C contains coefficients in the internal basis
2. **Householder Transformations**: Apply V_i = I - 2C_iC_i† to each photon i
3. **Mode Expansion**: Work in expanded space of n×m modes (n photons × m physical modes)
4. **Mode Shuffling**: Reorder modes from photon-major to mode-major ordering
5. **Physical Interferometer**: Apply m×m interferometer U to each photon's modes
6. **Clifford Sampling**: Use standard Clifford algorithm in expanded bosonic space
7. **Binning**: Collapse back to m physical modes by summing over photon blocks

## Files Created

### Core Implementation
- `src/boson_samplers/partial_distinguishability_householder.jl`
  Main sampler implementation with all helper functions

### Integration
- Modified `src/BosonSampling.jl`
  Added include for new sampler file

- Modified `src/boson_samplers/sample.jl`
  Added dispatch case for `UserDefinedGramMatrix` input type

### Examples and Tests
- `test/test_householder_integration.jl`
  Integration test verifying the sampler works with `sample!()`

- `examples/householder_sampler_usage.jl`
  Usage examples showing different scenarios

### Prototype
- `src/boson_samplers/partial_distinguishability_exact.jl`
  Original working prototype (not integrated, for reference)

## Usage

```julia
using BosonSampling

# Define system
n = 3  # photons
m = 5  # modes

# Create Gram matrix (positive semi-definite, Hermitian)
S = rand_gram_matrix_from_orthonormal_basis(n, 2)  # rank 2

# Create input with partial distinguishability
input = Input{UserDefinedGramMatrix}(first_modes(n, m), S)

# Create interferometer and event
interf = RandHaar(m)
ev = Event(input, FockSample(), interf)

# Sample using integrated interface
sample!(ev)

# Get result
output = ev.output_measurement.s
println("Sampled output: ", output)
```

## Key Functions

### `householder_sampler(ev::Event{UserDefinedGramMatrix, FockSample})`
Main entry point, called automatically by `sample!()` for `UserDefinedGramMatrix` inputs.

### `build_householder_interferometer(C, U, n, m, r)`
Constructs the full transformation matrix combining:
- Splitting matrices (Householder transformations)
- Shuffle permutation
- Block diagonal interferometers
- Reverse shuffle

### `build_splitting_matrices(C, n, m, r)`
Creates block diagonal matrix with Householder transformation for each photon.

### `build_shuffle_permutation(n, m)`
Creates permutation matrix for mode reordering (photon-major → mode-major).

### `build_block_diagonal_interferometers(U, n, m)`
Creates interleaved block diagonal structure with n copies of U.

### `bin_to_physical_modes(sampled_expanded, n, m)`
Collapses expanded n×m modes back to m physical modes.

## Gram Matrix Requirements

The Gram matrix S must satisfy:
- **Hermitian**: S = S†
- **Positive semi-definite**: all eigenvalues ≥ 0
- **Normalized**: S[i,i] = 1 for all photons i
- **Physical**: S[i,j] = ⟨photon_i, photon_j⟩ represents overlap

### Interpretation
- S[i,j] = 1: photons i and j are identical
- S[i,j] = 0: photons i and j are fully distinguishable
- 0 < |S[i,j]| < 1: partial distinguishability
- rank(S) = dimension of internal Hilbert space

## Advantages

1. **General**: Works with arbitrary Gram matrices (not limited to parameterized models)
2. **Exact**: No approximations in the sampling (uses Clifford algorithm exactly)
3. **Flexible**: Supports heterogeneous distinguishability between different photon pairs
4. **Integrated**: Works seamlessly with existing `sample!()` interface
5. **Compatible**: Fully compatible with existing Event/Input/Interferometer types

## Performance

Complexity: O(n²m² + Clifford_sampling)
- Scales with expanded mode space n×m
- Householder matrices: O(n×m²)
- Shuffle permutation: O((nm)²)
- Clifford sampling: dominant for large systems

## Validation

The implementation has been tested and verified to:
- Correctly preserve photon number (sum of output = n)
- Work with different Gram matrix ranks
- Integrate properly with `sample!()` dispatch
- Handle edge cases (fully indistinguishable, fully distinguishable)

## Future Enhancements

Potential improvements:
- Sparse matrix optimization for large m
- Precomputation and caching of transformation matrices
- Support for loss and detector imperfections
- Probability calculation via permanent computation

## References

- Householder transformation: https://en.wikipedia.org/wiki/Householder_transformation
- Gram matrix decomposition for partial distinguishability
- Clifford's corrected algorithm for boson sampling
