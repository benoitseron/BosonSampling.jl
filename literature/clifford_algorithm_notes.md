# Clifford & Clifford Algorithm Notes

From arXiv:1706.01260v2

## IMPORTANT: Matrix Convention in BosonSampling.jl

**This package uses a NON-STANDARD matrix convention that differs from the Clifford paper.**

### Standard Physics Convention (used in paper)
- U[i,j] = amplitude from input mode j to output mode i
- Rows = output modes, Columns = input modes
- Scattering matrix: S = U[output_modes, input_modes]

### BosonSampling.jl Convention
- U[i,j] = amplitude from input mode i to output mode j (indices swapped!)
- Rows = input modes, Columns = output modes
- Scattering matrix: S = U[input_modes, output_modes]

### Why This Matters
- For collision-free cases (each mode has ≤1 photon), both conventions give identical probabilities
  because the selected submatrices are transposes and Per(A) = Per(A^T)
- For cases with collisions (bunching), the conventions give DIFFERENT probabilities
- The Clifford sampler must transpose U before running to match compute_probability!

### Fix Applied
In `cliffords_sampler()`, we use:
```julia
A = transpose(interf.U)[:, fill_arrangement(input)]
```
instead of:
```julia
A = interf.U[:, fill_arrangement(input)]
```

This ensures the sampler produces distributions consistent with all other package functions.

## Problem Setup

- A = m×n complex matrix (first n columns of m×m Haar random unitary)
- m = output modes (rows)
- n = input photons (columns)
- Output z = (z₁,...,zₙ) sorted multiset

## Probability Formula

For sorted output z:
```
q(z) = |Per A_z|² / μ(z)
```
where:
- A_z is n×n matrix with rows z₁,...,zₙ from A
- μ(z) = ∏ sⱼ! (product of factorials of multiplicities)

## Expanded Sample Space

Sample r = (r₁,...,rₙ) from:
```
p(r) = (1/n!) |Per A_r|²
```
Then sort r to get z.

## Algorithm B (Fast O(n2ⁿ + poly(m,n)))

```
1: r ← ∅
2: A ← PERMUTE(A)                    # Randomly permute columns
3: wᵢ ← |aᵢ,₁|², i ∈ [m]
4: x ← SAMPLE(w)
5: r ← (r, x)
6: FOR k ← 2 TO n DO
7:    B_k^⋄ ← A_r^[k]                # (k-1)×k matrix: rows r, cols 1:k
8:    COMPUTE {Per B_k,ℓ^⋄, ℓ ∈ [k]} # minors: perm of (k-1)×(k-1) submatrices
9:    wᵢ ← |∑_{ℓ=1}^k aᵢ,ℓ Per B_k,ℓ^⋄|², i ∈ [m]  # Laplace expansion
10:   x ← SAMPLE(w)
11:   r ← (r, x)
12: END FOR
13: z ← INCSORT(r)
14: RETURN z
```

## Key Points

1. **Column permutation**: Essential for unbiased sampling. Random α allows averaging over permutations.

2. **Minors computation**: At step k, B_k^⋄ = A[r, 1:k] is (k-1)×k matrix.
   - Per B_k,ℓ^⋄ = permanent of B_k with column ℓ removed
   - This is a (k-1)×(k-1) square matrix

3. **Laplace expansion**: Weight for mode i at step k is:
   ```
   wᵢ = |∑_{ℓ=1}^k A[i,ℓ] × Per(B_k with col ℓ removed)|²
   ```
   This equals |Per(A[r∪{i}, 1:k])|² using Laplace expansion along the new row.

4. **Remark**: For Haar random unitary, columns already random, so can use identity permutation for single sample. BUT for multiple samples from same unitary, must re-permute each time.

## Matrix Convention

- Paper: A[i,j] = amplitude from input j to output i
- Rows index output modes
- Columns index input photons
