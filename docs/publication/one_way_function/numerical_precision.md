# Numerical Precision Analysis

> Precision model for the one-way-function CDF estimator implemented in
> [`src/boson_samplers/one_way_function.jl`](../../../src/boson_samplers/one_way_function.jl).
> The `F64-A` … `F64-E` tags in that source refer to the categories analyzed here.

## Overview

The CDF estimator $\hat{S}(x_0) = \frac{1}{M}\sum_{i=1}^M Z_i$ combines exact integer arithmetic (modular products via `BigInt`/`Int128`) with `Float64` complex arithmetic (phases, matrix products, geometric sums). We analyze whether the floating-point operations preserve the theoretical $O(\log N / \sqrt{M})$ convergence rate.

Each Monte Carlo sample:

1. Draws $k \sim q(k) \propto 1/\min(k, N-k)$ via rejection sampling on a continuous Float64 proposal
2. Draws $\mathbf{x} \sim \mathrm{Uniform}(\{-1,1\}^n)$
3. Computes $Z = \mathrm{Gly}_{\mathbf{x}}(U_{\mathrm{in}}^\dagger D_k U_{\mathrm{in}}) \cdot G_N(k, x_0) \;/\; (N \cdot q(k))$

We analyze each Float64 operation in turn and conclude with the one non-trivial concern: the Fourier mode sampling for $N > 2^{53}$.

---

## 1. Phase Computation

`_cispi2(r, N)` computes $e^{2\pi i r/N}$ via `cispi(2 \cdot \texttt{Float64}(r) / \texttt{Float64}(N))$.

**Supplementary angle trick.** When $r > N/2$, the code evaluates $e^{-2\pi i (N-r)/N}$ instead. This ensures the argument to `cispi` satisfies $|x| \leq 1$, preventing the degeneracy where $\texttt{Float64}(r)/\texttt{Float64}(N)$ rounds to $1.0$ for $r \approx N$.

**Error bound.** The ratio $r/N \in [0, 1/2]$ has absolute Float64 error $\leq 2\varepsilon_\mathrm{mach}$. Julia's `cispi` uses `sinpi`/`cospi` with exact argument reduction for half-integers, giving:

$$|\Delta e^{2\pi i r/N}| \leq 2\pi \cdot 2\varepsilon_\mathrm{mach} \approx 1.4 \times 10^{-15}$$

**Empirical validation.** Comparing against 256-bit BigFloat reference for $N$ up to $10^{65}$, the measured error is $\leq 3 \times 10^{-16}$ in all cases, including BigInt arguments.

---

## 2. Geometric Sum

For $k \neq 0$:

$$G_N(k, x_0) = \frac{1 - e^{2\pi i k(x_0+1)/N}}{1 - e^{2\pi i k/N}}$$

The phase indices $k(x_0+1) \bmod N$ and $k \bmod N$ are computed in exact modular arithmetic (`_mulmod`). Only the final phase evaluation uses Float64.

**Near-cancellation regime** ($k/N \ll 1$): the subtraction $1 - e^{i\theta}$ for small $\theta$ is dominated by its imaginary part $-\sin\theta \approx -\theta$, which `sinpi` computes to full precision. The real part $1 - \cos\theta$ loses bits but is subdominant. The ratio of two such terms is well-conditioned.

**Empirical validation.** Relative error $\leq 2.2 \times 10^{-16}$ ($\approx \varepsilon_\mathrm{mach}$) for $N$ up to $10^{65}$ (compared against BigFloat closed-form).

---

## 3. Fused Glynn Estimator

The fused evaluator computes $\mathrm{Gly}_\mathbf{x}(U_\mathrm{in}^\dagger D_k U_\mathrm{in})$ in $O(mn)$ operations without materializing the $n \times n$ matrix. Standard forward error analysis for the product of $n$ inner products, each involving $m$ terms:

$$|\Delta\mathrm{Gly}| \lesssim (2m + n)\,\varepsilon_\mathrm{mach}$$

since $|\mathrm{Gly}_\mathbf{x}(A)| \leq \|A\|^n \leq 1$ for subunitary $A$. For $(n,m) = (40,40)$: $|\Delta\mathrm{Gly}| \lesssim 1.2 \times 10^{-14}$.

**Empirical validation.** Over 200 random instances with $n \in [2,8]$, the fused evaluator matches the explicit-matrix Glynn estimator to within $3.7\,\varepsilon_\mathrm{mach}$.

---

## 4. Monte Carlo Accumulation

The running sum of $M$ complex samples, each bounded by $|Z| = O(\log^2 N)$ (Proposition 3), accumulates rounding error. After averaging:

$$|\Delta\hat{S}_\mathrm{round}| \lesssim \varepsilon_\mathrm{mach} \cdot \log^2 N$$

For $N \sim 10^{65}$: $\lesssim 10^{-16} \cdot 150^2 \approx 2 \times 10^{-12}$, negligible compared to the MC statistical error $O(\log N / \sqrt{M}) \sim 10^{-2}$.

---

## 5. Fourier Mode Sampling (Main Concern)

The importance sampler draws $k$ from $q(k) \propto 1/\min(k, N\!-\!k)$ by:

1. Choosing branch: $k = 0$ (mass $\propto 1$), lower ($k \in [1, N/2\!-\!1]$, mass $\propto H_{N/2-1}$), or upper ($k \in [N/2, N\!-\!1]$, mass $\propto H_{N-N/2}$).
2. Within each branch: rejection sampling on the Float64 proposal $k = \lfloor (K+1)^V \rfloor$, $V \sim \mathrm{Uniform}[0,1]$.

**The discretization issue.** Float64 has a 53-bit mantissa. A Float64 proposal can represent at most $\sim 2^{53}$ distinct values of $k$. All Fourier modes with $\min(k, N\!-\!k) \leq 2^{52}$ are reachable; beyond this threshold the proposal reaches only a sparse, $\sim 2^{53}$-point lattice of modes, so the **actual** sampling density $q'(k)$ differs from the **intended** $q(k) \propto 1/\min(k,N\!-\!k)$.

**This is a real, deterministic bias — not negligible (corrected).** Earlier versions of this section argued the effect away by bounding the *contribution of each dropped high-$k$ mode*. That reasoning is wrong: the estimator divides by the *intended* weight $1/q(k)$, so a sampler that draws from $q'(k) \neq q(k)$ gives

$$ \mathbb{E}[Z] = \sum_k \frac{q'(k)}{q(k)}\,\frac{\varphi(k)\,G_N(k,x_0)}{N} \;\neq\; S(x_0). $$

The mismatch *reweights the whole spectrum*, including the $O(1)$ low-$k$ contributions — it is not the sum of tiny dropped terms. Measured against base-blind ground truth (direct boson sampling; the $f$-ordering, hence the true $S$, is base-independent for any base $b \geq n+1$):

| case | true $S$ | $\widehat S$ (Float64 proposal) |
|------|----------|----------------------------------|
| $n{=}12, m{=}144$, base $17$ | $0.515$ | $0.54$ (ok) |
| $n{=}12, m{=}144$, base $16{=}2^4$ | $0.515$ | $\mathbf{1.16}$ |
| $n{=}12, m{=}144$, base $32{=}2^5$ | $0.515$ | $\mathbf{1.24}$ |
| $n{=}4, m{=}40$, base $16$ | $0.428$ | $\mathbf{0.83}$ |

The bias grows with $\mathrm{bits}(N) = m\log_2 b$ and is **catastrophic for power-of-two encoding bases**: then $N = n\,2^{a(m-1)}{+}1$ makes $K=N/2$ a few-significant-bit (dyadic) number, so the Float64 lattice is *structurally aligned* with the (also dyadic) spectral mass and the error is coherent rather than averaging out. Generic bases show the same effect far more slowly. The old "argmax stable / sum-to-one" checks missed it because they never compared the *absolute* $\widehat S$ at large $N$ to ground truth.

### Fix (implemented)

The proposal $k=\lfloor (K{+}1)^v\rfloor$ is now evaluated with $\gtrsim \mathrm{bits}(K)$ precision (`BigFloat` via `setprecision`) whenever $K > 2^{52}$, in **both** `_sample_reciprocal` (Int64/Int128) and `_sample_reciprocal!` (BigInt). At that precision $\lfloor\cdot\rfloor$ resolves the exact integer mode, so $q'(k)=q(k)$ and the estimator is unbiased again. The fast Float64 path is retained for $K \leq 2^{52}$ (already exact). Cost is one high-precision `pow` per sample, comparable to the per-sample BigInt arithmetic already performed.

**Verification.** With the fix, the four cases above return $\widehat S = 0.49, 0.49, 0.53, 0.43$ (all within MC noise of the truth), and the full test suite (2272 tests) passes. The earlier-suspected "power-of-two resonance" was entirely this sampler bug.

### Formal status

For $N \leq 2^{53}$: provably unbiased (all modes reachable, Float64 proposal exact). For $N > 2^{53}$: unbiased to the working precision $\gtrsim\mathrm{bits}(K)$ of the BigFloat proposal (residual $\lesssim 2^{-16}$ probability of a $\pm1$-mode misresolution, statistically harmless).

---

## 6. Note on Type Stability

`_compute_N` returns `Int64`, `Int128`, or `BigInt` depending on magnitude, which is inherently type-unstable. To prevent this union type from infecting the hot path, `find_most_probable_bin` uses a **function barrier**: the type-unstable `_compute_N` call is isolated in the outer function, while the inner loop (`_find_most_probable_bin`) receives the concrete type `T` via `where T <: Integer` and is fully specialized by the compiler.

Similarly, `estimate_S` dispatches to the zero-allocation BigInt path via method dispatch rather than a runtime `isa` check, and the hot-path functions (`_Z_sample!`, `_glynn_fused!`) are parametric on `T`.

Verified with `@code_warntype`: the inner function shows no `Union` type instabilities for any concrete integer type.

---

## Summary

| Operation | Per-call error | Regime | vs. MC error |
|-----------|---------------|--------|-------------|
| Phase $e^{2\pi i r/N}$ | $\leq 1.4 \times 10^{-15}$ | All $N$ | negligible |
| Geometric sum $G_N$ | $\leq \varepsilon_\mathrm{mach}$ (relative) | All $N$ | negligible |
| Glynn product | $\leq (2m+n)\varepsilon_\mathrm{mach}$ | All $N$ | negligible |
| MC accumulation | $\leq \varepsilon_\mathrm{mach} \log^2 N$ | All $N$ | negligible |
| Fourier sampling | **exact** for $K \leq 2^{52}$ (Float64) | $n+m \lesssim 35$ | zero |
| Fourier sampling | **exact** for $K > 2^{52}$ (BigFloat, $\gtrsim$bits$(K)$) | BigInt regime | negligible ($\lesssim 2^{-16}$) |

**Conclusion.** The Fourier-mode proposal is evaluated at precision sufficient to resolve the exact integer mode in every regime (Float64 for $K\leq2^{52}$, BigFloat above), so the sampler is unbiased and the only remaining error is the Monte Carlo variance $O(\log^2 N / M)$. *Historical note:* before this fix the Float64 proposal biased $\widehat S(x_0)$ at large $N$ — severely for power-of-two encoding bases — which earlier surfaced as a spurious "resonance" spike in the sample-budget figure.
