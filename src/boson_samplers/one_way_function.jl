# ═══════════════════════════════════════════════════════════════════════
# One-way function via boson sampling: CDF estimation & most-probable bin
# ═══════════════════════════════════════════════════════════════════════
# Integrated into BosonSampling as a module include(); dependencies
# (LinearAlgebra, Random, StatsBase, Permanents, and RandHaar) are provided at
# module scope by src/BosonSampling.jl.  Stars-and-bars counts use Base
# `binomial(big(...), big(...))` (exact BigInt) — deliberately NOT
# BigCombinatorics.Binomial, whose name would shadow Distributions.Binomial.
#
# Reproducible paper figures and speed tests live in
# docs/publication/one_way_function/ (see numerical_precision.md for the
# floating-point precision model referenced by the F64-<X> tags below).
#
# No `export` here: BosonSampling re-exports every module identifier via the
# auto-export loop at the end of src/BosonSampling.jl, so `using BosonSampling`
# already brings the whole API (and internals) into scope.

# ═══════════════════════════════════════════════════════════════════════
# Floating-point precision model
# ═══════════════════════════════════════════════════════════════════════
# Integer arithmetic (phase indices, modular products) is exact via
# Int64 / Int128 / BigInt.  Five categories of operations approximate
# via Float64.  Each is marked inline with "# F64-<X>" for traceability.
#
# F64-A  Phase computation: cispi(2·Float64(r)/Float64(N)).
#        Supplementary angle trick keeps the argument in [-1, 1].
#        Per-call error: |Δe^{2πir/N}| ≤ 2π·2ε_mach ≈ 1.4×10⁻¹⁵.
#        Exact when N ≤ 2⁵³; ε_mach relative error otherwise.
#
# F64-B  Importance weight: t(k) = 1/Float64(k).
#        Relative error ε_mach.  The same Float64(k) value appears in
#        both q(k) and the 1/q(k) correction, so the errors cancel
#        exactly in the estimator Z = Gly·G_N/(N·q(k)).
#        No bias contribution.
#
# F64-C  Fourier-mode proposal: k = floor((K+1)^v).  A Float64 proxy resolves
#        only ~53 bits, so for K>2⁵² it reaches a sparse lattice of modes and
#        samples the WRONG discrete distribution q'(k)≠q(k) — a real, deterministic
#        BIAS on Ŝ(x₀) (not negligible; catastrophic for power-of-two encoding
#        bases).  FIXED: the proposal is now evaluated with ≳bits(K) precision
#        (BigFloat) for K>2⁵², which reaches the exact mode and removes the bias.
#        See docs/publication/one_way_function/numerical_precision.md §5.
#
# F64-D  Normalization divisor: Z is divided by Float64(N).
#        Introduces multiplicative bias (1 ± ε_mach) on Ŝ(x₀), which
#        cancels in bin-probability differences p̂_j = Ŝ(x_{j+1})−Ŝ(x_j).
#
# F64-E  Type-selection guard in _compute_N: Float64 magnitude comparison
#        only; not used in any downstream computation.
#
# With the F64-C proposal evaluated at adequate precision (see below), all
# approximation errors are negligible vs. the MC statistical error O(log²N / M).
# See docs/publication/one_way_function/numerical_precision.md for details.
# ═══════════════════════════════════════════════════════════════════════

# ─── Glynn single-sample estimator ────────────────────────────────────
# Returns Gly_x(A) = prod(x) * prod_i(sum_j A[i,j]*x[j])
# for a single random x ∈ {-1,1}^n (Eq. 14 of the paper / Eq. 20 of Aaronson-Hance).
function glynn_single(A::AbstractMatrix, x::AbstractVector)
    n = size(A, 1)
    val = prod(x)
    for i in 1:n
        val *= sum(A[i, j] * x[j] for j in 1:n)
    end
    return val
end

# ─── Importance sampling distribution q(k) (Proposition 2) ───────────
# t(k) as in Eq. 17, with N Fourier modes indexed 0:N-1.
function t_weight(k::Integer, N::Integer)
    if k == 0
        return 1.0
    elseif 1 <= k <= N ÷ 2 - 1
        return 1.0 / Float64(k)       # F64-B: cancels with 1/q(k) in estimator
    else  # N/2 <= k <= N-1
        return 1.0 / Float64(N - k)   # F64-B
    end
end

# log(k) as Float64, accurate even when Float64(k) would overflow to Inf
# (BigInt N can exceed the Float64 ceiling ~1.8e308, i.e. n·(n+1)^(m-1) ≳ 1e308).
@inline function _log_int(k::Integer)
    kf = Float64(k)
    isfinite(kf) && return log(kf)
    # k too large for Float64: log(k) = e·log(2) + log(k >> e), keeping ~60 mantissa bits.
    e = ndigits(k, base = 2) - 60
    return e * 0.6931471805599453 + log(Float64(k >> e))
end

# Harmonic number H_k = sum_{i=1}^k 1/i.
# For large k, uses asymptotic expansion (error < 10^{-15} for k > 10000).
function _harmonic(k::Integer)
    if k <= 0
        return 0.0
    elseif k <= 10_000
        s = 0.0
        for i in 1:k
            s += 1.0 / Float64(i)  # Float64(i): exact for i ≤ 10000; type-stable for BigInt k
        end
        return s
    else
        γ = 0.5772156649015329
        L = _log_int(k)  # F64-B: H_k only enters via importance weights; overflow-safe
        x = Float64(k)
        isfinite(x) || return L + γ  # 1/x… correction terms < 1e-300 when k > 1.8e308
        x2 = x * x
        x4 = x2 * x2
        x6 = x4 * x2
        return L + γ + 1 / (2x) - 1 / (12x2) + 1 / (120x4) - 1 / (252x6)
    end
end

# Normalization constant 𝒩 = sum of t(k) for k=0..N-1  (Eq. 22).
# Uses harmonic numbers: 𝒩 = 1 + H_{⌊N/2⌋-1} + H_{N-⌊N/2⌋}.
function normalization_constant(N::Integer)
    half = N ÷ 2
    return 1.0 + _harmonic(half - 1) + _harmonic(N - half)
end

# ─── Precomputed sampling context ────────────────────────────────────
# Avoids recomputing O(N) harmonic sums on every MC sample.
# Parametric in the integer type T so that N can be Int64, Int128, or BigInt.
struct SamplingContext{T<:Integer}
    N::T
    half::T
    half_m1::T          # half - 1 (precomputed to avoid BigInt alloc)
    upper_range::T      # N - half (precomputed)
    mass_zero::Float64
    mass_lower::Float64
    mass_upper::Float64
    total::Float64      # = 𝒩
end

function SamplingContext(N::Integer)
    half = N ÷ 2
    half_m1 = half - 1
    upper_range = N - half
    mass_zero = 1.0
    mass_lower = _harmonic(half_m1)
    mass_upper = _harmonic(upper_range)
    total = mass_zero + mass_lower + mass_upper
    return SamplingContext(N, half, half_m1, upper_range, mass_zero, mass_lower, mass_upper, total)
end

# Sample k from q(k) = t(k)/𝒩 using precomputed context.
function sample_fourier_mode(ctx::SamplingContext)
    u = rand() * ctx.total
    if u < ctx.mass_zero
        return zero(ctx.N)
    elseif u < ctx.mass_zero + ctx.mass_lower
        return _sample_reciprocal(ctx.half_m1)
    else
        kprime = _sample_reciprocal(ctx.upper_range)
        return ctx.N - kprime
    end
end

# Legacy interface (recomputes masses — for testing only).
function sample_fourier_mode(N::Integer)
    return sample_fourier_mode(SamplingContext(N))
end

# Sample an integer k ∈ {1,...,K} with probability ∝ 1/k using the inverse-CDF
# proxy floor((K+1)^v), v~U[0,1] (CDF F_Y(k)=log(k+1)/log(K+1)), plus rejection
# (Eq. 18-21).  The proxy must be evaluated with enough precision to resolve a
# UNIQUE integer in [1,K].  A Float64 proxy (~53-bit mantissa) cannot do this for
# K>2⁵²: it reaches only a sparse, ~2⁵³-point lattice of k, sampling the wrong
# discrete distribution q'(k)≠q(k) and biasing the importance-weighted estimator
# (catastrophically for power-of-two encoding bases; see
# docs/publication/one_way_function/numerical_precision.md §5).  We keep the fast Float64 path for
# K≤2⁵² (exact, all integers reachable) and use a ≳bits(K)-precision BigFloat
# proxy above it.
function _sample_reciprocal(K::Integer)
    T = typeof(K)
    if K <= (one(K) << 52)          # Float64 exact: every integer in [1,K] reachable
        Kf = Float64(K)
        while true
            v = rand()
            k = clamp(floor(T, (Kf + 1.0)^v), oneunit(T), K)
            kf = Float64(k)
            rand() <= log(2) / (kf * log1p(1.0 / kf)) && return k
        end
    else                            # K>2⁵²: resolve the exact mode with full precision
        prec = ndigits(K, base = 2) + 16
        return setprecision(BigFloat, prec) do
            Kp1 = BigFloat(K) + 1
            while true
                v = rand()
                k = clamp(floor(T, Kp1^v), oneunit(T), K)
                kf = Float64(k)
                acc = isfinite(kf) ? log(2) / (kf * log1p(1.0 / kf)) : 0.6931471805599453
                rand() <= acc && return k
            end
        end
    end
end

# ─── Overflow-safe modular multiplication ────────────────────────────
# a * b mod N, using a wider type to prevent overflow.
@inline _mulmod(a::Int64, b::Int64, N::Int64) = Int64(mod(Int128(a) * Int128(b), Int128(N)))
@inline _mulmod(a::Int128, b::Int128, N::Int128) = Int128(mod(big(a) * big(b), big(N)))
@inline _mulmod(a::BigInt, b::BigInt, N::BigInt) = mod(a * b, N)

# In-place modular multiplication for BigInt: r = (r * b) mod N using preallocated temp.
# Avoids all allocation in the hot loop.
@inline function _mulmod!(r::BigInt, b::BigInt, N::BigInt, temp::BigInt)
    Base.GMP.MPZ.mul!(temp, r, b)
    Base.GMP.MPZ.tdiv_r!(r, temp, N)
    return r
end

# ─── Phase helper for large integers ─────────────────────────────────
# Compute cispi(2r/N) = exp(2πi r/N) accurately for r ∈ [0, N).
# When r is close to N, Float64(r)/Float64(N) rounds to 1.0 and cispi(2.0)=1
# gives a spurious zero. Fix: use the supplementary angle when r > N/2.
# half: precomputed N ÷ 2 to avoid BigInt allocation per call.
@inline function _cispi2(r::Integer, N::Integer, half::Integer)
    if r <= half          # F64-A: r/N ∈ [0,½], error ≤ 2π·2ε_mach per call
        return cispi(2.0 * Float64(r) / Float64(N))
    else                  # F64-A: supplementary angle avoids r/N ≈ 1 degeneracy
        return cispi(-2.0 * Float64(N - r) / Float64(N))
    end
end
@inline _cispi2(r::Integer, N::Integer) = _cispi2(r, N, N ÷ 2)

# ─── Geometric sum G_N(k; x0) (Eq. 9) ────────────────────────────────
function G_N(k::Integer, x0::Integer, N::Integer)
    if k == 0
        return Complex(Float64(x0) + 1.0)  # F64-D: exact for x0 < 2⁵³
    else
        # Use _mulmod to prevent overflow in k*(x0+1), then _cispi2 for phase precision.
        num_phase = _mulmod(mod(k, N), mod(x0 + 1, N), N)
        den_phase = mod(k, N)
        return (1 - _cispi2(num_phase, N)) / (1 - _cispi2(den_phase, N))
    end
end

# ─── Diagonal phase matrix D_k (for testing) ─────────────────────────
# D_k = diag(exp(-2πi k ω_1/N), ..., exp(-2πi k ω_m/N))
function diag_Dk(k::Integer, omega::AbstractVector, N::Integer)
    return Diagonal([conj(_cispi2(mod(big(k) * round(BigInt, ω_j), N), N)) for ω_j in omega])
end

# ─── Fused Glynn + phase evaluation ───────────────────────────────────
# Computes Gly_x(U_in† D_k U_in) in O(m·n) without materializing the n×n matrix.
# Uses modular arithmetic for phase indices to avoid overflow.
function _glynn_fused!(buf::AbstractVector, U_in::AbstractMatrix,
                       k::T, base::T, N::T, half::T,
                       x::Vector{Int}) where T <: Integer
    C = eltype(buf)
    m, n = size(U_in)

    @inbounds begin
        xj = C(x[1])
        for l in 1:m; buf[l] = U_in[l, 1] * xj; end
        for j in 2:n
            xj = C(x[j])
            for l in 1:m; buf[l] += U_in[l, j] * xj; end
        end
    end

    r = mod(k, N)
    @inbounds for l in 1:m
        buf[l] *= conj(_cispi2(r, N, half))
        r = _mulmod(r, base, N)
    end

    result = C(prod(x))
    @inbounds for i in 1:n
        row_sum = zero(C)
        for l in 1:m; row_sum += conj(U_in[l, i]) * buf[l]; end
        result *= row_sum
    end
    return result
end

# ─── Single MC sample ────────────────────────────────────────────────
function _Z_sample!(buf::AbstractVector, x::Vector{Int},
                    U_in::AbstractMatrix, base::T, N::T,
                    x0::T, ctx::SamplingContext{T}) where T <: Integer
    k = sample_fourier_mode(ctx)
    @inbounds for j in eachindex(x)
        x[j] = ifelse(rand(Bool), 1, -1)
    end
    gly = _glynn_fused!(buf, U_in, k, base, N, ctx.half, x)
    g = G_N(k, x0, N)
    qk = t_weight(k, N) / ctx.total
    return gly * g / (Float64(N) * qk)  # F64-D: bias (1±ε_mach), cancels in differences
end

# Legacy Z_sample (for testing — not optimized).
function Z_sample(U::AbstractMatrix, omega::AbstractVector, N_fourier::Integer,
                  x0::Integer, n_photons::Int)
    k = sample_fourier_mode(N_fourier)
    Dk = diag_Dk(k, omega, N_fourier)
    A = U' * Dk * U
    A_sub = A[1:n_photons, 1:n_photons]
    x = rand([-1, 1], n_photons)
    gly = glynn_single(A_sub, x)
    𝒩 = normalization_constant(N_fourier)
    qk = t_weight(k, N_fourier) / 𝒩
    return gly * G_N(k, x0, N_fourier) / (N_fourier * qk)
end

# ═══════════════════════════════════════════════════════════════════════
# Zero-allocation BigInt hot path
# ═══════════════════════════════════════════════════════════════════════
# When N is BigInt, every arithmetic op allocates GMP integers.
# This path uses preallocated workspace and in-place GMP ops (MPZ.mul!, etc.)
# to eliminate per-sample allocations, enabling linear thread scaling.

struct _BigIntWork
    k::BigInt       # sampled Fourier mode
    kprime::BigInt  # scratch for upper branch of sample_fourier_mode
    r::BigInt       # phase recurrence in Glynn loop
    temp::BigInt    # scratch for _mulmod! and N-r subtraction
    temp2::BigInt   # scratch for G_N numerator phase
end
_BigIntWork() = _BigIntWork(BigInt(), BigInt(), BigInt(), BigInt(), BigInt())

# floor(Float64) → BigInt, in-place. GMP's mpz_set_d truncates toward zero,
# which equals floor for positive values (always the case here).
@inline function _set_floor_pos!(z::BigInt, x::Float64)
    ccall((:__gmpz_set_d, :libgmp), Cvoid, (Base.GMP.MPZ.mpz_t, Cdouble), z, x)
    return z
end

# r/N ∈ [0,1) as Float64, accurate to ~53 bits, robust when Float64(N) overflows.
# r assumed in [0, N).  Used for phase ratios and importance weights in the BigInt path.
@inline function _ratio(r::BigInt, N::BigInt, Nf::Float64)
    isfinite(Nf) && return Float64(r) / Nf
    return Float64((r << 53) ÷ N) / 9.007199254740992e15  # (r·2^53 ÷ N) / 2^53
end

# Sample k ∈ {1,...,K} with prob ∝ 1/k, writing result into ws.k. Zero-alloc when
# K fits Float64; falls back to log-space sampling (allocating) for huge BigInt K.
function _sample_reciprocal!(out::BigInt, K::BigInt)
    if K <= (BigInt(1) << 52)       # Float64 exact: fast path, all integers reachable
        Kf = Float64(K)
        while true
            v = rand()
            _set_floor_pos!(out, (Kf + 1.0)^v)
            out < 1   && Base.GMP.MPZ.set_si!(out, 1)
            out > K   && Base.GMP.MPZ.set!(out, K)
            kf = Float64(out)
            rand() <= log(2) / (kf * log1p(1.0 / kf)) && return out
        end
    else                            # K>2⁵²: ≳bits(K)-precision proxy resolves the exact mode
        prec = ndigits(K, base = 2) + 16
        return setprecision(BigFloat, prec) do
            Kp1 = BigFloat(K) + 1
            while true
                v = rand()
                Base.GMP.MPZ.set!(out, floor(BigInt, Kp1^v))
                out < 1   && Base.GMP.MPZ.set_si!(out, 1)
                out > K   && Base.GMP.MPZ.set!(out, K)
                kf = Float64(out)
                acc = isfinite(kf) ? log(2) / (kf * log1p(1.0 / kf)) : 0.6931471805599453
                rand() <= acc && return out
            end
        end
    end
end

# Sample Fourier mode into ws.k. Zero-alloc.
function _sample_fourier_mode!(ws::_BigIntWork, ctx::SamplingContext{BigInt})
    u = rand() * ctx.total
    if u < ctx.mass_zero
        Base.GMP.MPZ.set_si!(ws.k, 0)
    elseif u < ctx.mass_zero + ctx.mass_lower
        _sample_reciprocal!(ws.k, ctx.half_m1)
    else
        _sample_reciprocal!(ws.kprime, ctx.upper_range)
        Base.GMP.MPZ.sub!(ws.k, ctx.N, ws.kprime)
    end
    return ws.k
end

# Fused Glynn evaluator for BigInt N. Zero per-sample allocation.
function _glynn_fused_bigint!(buf::AbstractVector, U_in::AbstractMatrix,
                               ws::_BigIntWork, base::BigInt, N::BigInt,
                               half::BigInt, Nf::Float64, x::Vector{Int})
    C = eltype(buf)
    m, n = size(U_in)

    @inbounds begin
        xj = C(x[1])
        for l in 1:m; buf[l] = U_in[l, 1] * xj; end
        for j in 2:n
            xj = C(x[j])
            for l in 1:m; buf[l] += U_in[l, j] * xj; end
        end
    end

    # Phase recurrence: r = k, k·base, k·base², ... (mod N)
    Base.GMP.MPZ.set!(ws.r, ws.k)
    # F64-A: phase recurrence uses Float64 ratio r/N (supplementary angle trick)
    @inbounds for l in 1:m
        if ws.r <= half
            buf[l] *= cispi(-2.0 * _ratio(ws.r, N, Nf))
        else
            Base.GMP.MPZ.sub!(ws.temp, N, ws.r)
            buf[l] *= cispi(2.0 * _ratio(ws.temp, N, Nf))
        end
        _mulmod!(ws.r, base, N, ws.temp)
    end

    result = C(prod(x))
    @inbounds for i in 1:n
        row_sum = zero(C)
        for l in 1:m; row_sum += conj(U_in[l, i]) * buf[l]; end
        result *= row_sum
    end
    return result
end

# Single MC sample for BigInt N. Zero per-sample allocation.
# x0p1_mod_N = mod(x0 + 1, N), precomputed once per estimate_S call.
function _Z_sample_bigint!(buf::AbstractVector, x::Vector{Int},
                            U_in::AbstractMatrix, base::BigInt, N::BigInt,
                            half::BigInt, Nf::Float64, x0f::Float64,
                            x0p1_mod_N::BigInt,
                            ctx::SamplingContext{BigInt}, ws::_BigIntWork)
    _sample_fourier_mode!(ws, ctx)

    @inbounds for j in eachindex(x)
        x[j] = ifelse(rand(Bool), 1, -1)
    end

    gly = _glynn_fused_bigint!(buf, U_in, ws, base, N, half, Nf, x)

    # Term = gly·G_N·𝒩·ρ, where ρ = t_weight·(distance to nearer endpoint)/N ∈ (0,½].
    # Writing every "1/N" as a bounded ratio (via _ratio) keeps quantities finite when N
    # exceeds the Float64 ceiling.  For finite Nf this is algebraically identical to the
    # original gly·G_N/(N·qk).
    if iszero(ws.k)
        # G_N(0) = x0+1, t_weight(0) = 1  ⇒  term = gly·𝒩·((x0+1)/N).
        # x0p1_mod_N == 0 means x0+1 = N exactly (x0 = N-1), i.e. ratio = 1.
        r0 = iszero(x0p1_mod_N) ? 1.0 : _ratio(x0p1_mod_N, N, Nf)
        return gly * ctx.total * r0
    end

    # G_N numerator phase: A = 2π·(k(x0+1) mod N)/N (general angle, well-conditioned).
    Base.GMP.MPZ.set!(ws.temp2, x0p1_mod_N)
    _mulmod!(ws.temp2, ws.k, N, ws.temp)
    if ws.temp2 <= half
        num_cis = cispi(2.0 * _ratio(ws.temp2, N, Nf))
    else
        Base.GMP.MPZ.sub!(ws.temp, N, ws.temp2)
        num_cis = cispi(-2.0 * _ratio(ws.temp, N, Nf))
    end

    # Combine the denominator 1-e^{iB} (B = 2π·k/N) with the weight ρ = d/N as a single
    # factor W = ρ/(1-e^{iB}).  The importance sampler concentrates on d = min(k,N-k) small,
    # so when N is astronomically large B underflows and 1-e^{iB} rounds to 0 — but the limit
    # ρ/(1-e^{iB}) → ±i/(2π) is finite (the d/N factors cancel).  (F64-A: supplementary angle.)
    if ws.k <= half
        ρ = _ratio(ws.k, N, Nf)                       # d = k
        denom = 1 - cispi(2.0 * ρ)
        W = iszero(denom) ? Complex(0.0, 0.15915494309189535) : ρ / denom   # +i/(2π)
    else
        Base.GMP.MPZ.sub!(ws.temp, N, ws.k)           # d = N-k
        ρ = _ratio(ws.temp, N, Nf)
        denom = 1 - cispi(-2.0 * ρ)
        W = iszero(denom) ? Complex(0.0, -0.15915494309189535) : ρ / denom  # -i/(2π)
    end
    return gly * (1 - num_cis) * ctx.total * W
end

# ─── Estimate cumulative distribution S(x0) ──────────────────────────
# Generic path for Int64 / Int128.  BigInt dispatches to zero-alloc specialization below.
function estimate_S(U_in::AbstractMatrix, base::T, N::T,
                    x0::T, n::Int, M::Int, ctx::SamplingContext{T}) where T <: Integer
    C = eltype(U_in)
    m = size(U_in, 1)
    nt = Threads.nthreads()
    chunk = cld(M, nt)
    partials = Vector{C}(undef, nt)
    Threads.@threads for tid in 1:nt
        buf = Vector{C}(undef, m)
        x = Vector{Int}(undef, n)
        local_total = zero(C)
        local_M = min(chunk, M - (tid - 1) * chunk)
        for _ in 1:local_M
            local_total += _Z_sample!(buf, x, U_in, base, N, x0, ctx)
        end
        partials[tid] = local_total
    end
    return real(sum(partials) / M)
end

# BigInt specialization: zero per-sample allocation via preallocated workspace.
function estimate_S(U_in::AbstractMatrix, base::BigInt, N::BigInt,
                    x0::BigInt, n::Int, M::Int, ctx::SamplingContext{BigInt})
    _estimate_S_bigint(U_in, base, N, x0, n, M, ctx)
end

function _estimate_S_bigint(U_in::AbstractMatrix, base::BigInt, N::BigInt,
                             x0::BigInt, n::Int, M::Int, ctx::SamplingContext{BigInt})
    C = eltype(U_in)
    m = size(U_in, 1)
    Nf = Float64(N)    # F64-A,D: precomputed; used for phases and normalization
    x0f = Float64(x0)  # F64-D: only used in G_N(0) = x0+1 shortcut
    half = ctx.half
    x0p1_mod_N = mod(x0 + 1, N)  # precompute once (not per-sample)
    nt = Threads.nthreads()
    chunk = cld(M, nt)
    partials = Vector{C}(undef, nt)
    Threads.@threads for tid in 1:nt
        buf = Vector{C}(undef, m)
        x = Vector{Int}(undef, n)
        ws = _BigIntWork()  # one workspace per thread
        local_total = zero(C)
        local_M = min(chunk, M - (tid - 1) * chunk)
        for _ in 1:local_M
            local_total += _Z_sample_bigint!(buf, x, U_in, base, N, half, Nf,
                                              x0f, x0p1_mod_N, ctx, ws)
        end
        partials[tid] = local_total
    end
    return real(sum(partials) / M)
end

# Raw single-shot samples: returns the vector of K individual real Z-estimates that
# estimate_S would otherwise average away. Use this to estimate the per-shot variance
# σ² = Var(one Z-sample) directly — one sample-variance over K shots is far cheaper than
# nesting an inner mean inside outer repeats (see paper_figure.jl panel (d)).
function z_samples(U_in::AbstractMatrix, base::T, N::T,
                   x0::T, n::Int, K::Int, ctx::SamplingContext{T}) where T <: Integer
    m = size(U_in, 1)
    out = Vector{Float64}(undef, K)
    nt = Threads.nthreads()
    chunk = cld(K, nt)
    Threads.@threads for tid in 1:nt
        buf = Vector{eltype(U_in)}(undef, m)
        x = Vector{Int}(undef, n)
        for i in ((tid - 1) * chunk + 1):min(tid * chunk, K)
            out[i] = real(_Z_sample!(buf, x, U_in, base, N, x0, ctx))
        end
    end
    return out
end

# BigInt specialization (mirrors _estimate_S_bigint's per-thread workspace).
function z_samples(U_in::AbstractMatrix, base::BigInt, N::BigInt,
                   x0::BigInt, n::Int, K::Int, ctx::SamplingContext{BigInt})
    m = size(U_in, 1)
    Nf = Float64(N); x0f = Float64(x0); half = ctx.half
    x0p1_mod_N = mod(x0 + 1, N)
    out = Vector{Float64}(undef, K)
    nt = Threads.nthreads()
    chunk = cld(K, nt)
    Threads.@threads for tid in 1:nt
        buf = Vector{eltype(U_in)}(undef, m)
        x = Vector{Int}(undef, n)
        ws = _BigIntWork()
        for i in ((tid - 1) * chunk + 1):min(tid * chunk, K)
            out[i] = real(_Z_sample_bigint!(buf, x, U_in, base, N, half, Nf,
                                            x0f, x0p1_mod_N, ctx, ws))
        end
    end
    return out
end

# ────────────────────────────────────────────────────────────────────
# Fast Clifford & Clifford 2018 boson sampler
# ────────────────────────────────────────────────────────────────────
# Reimplemented to avoid the perf pitfalls of BosonSampling.cliffords_sampler:
# no `global` state, no `Threads.@threads` launch on tiny inner permanents
# (thread overhead was dominating for small n), pre-allocated workspace, and
# single-pass cumulative-weight sampling.  Convention identical to
# `exact_probability`: Pr[s] = |Per(U[s_modes, 1:n])|² / ∏ s_j! for input
# |1ⁿ 0^(m−n)⟩ — verified by 50k-shot empirical test against the exact PMF.
struct CCSamplerWorkspace
    A::Matrix{ComplexF64}            # m × n  (= U[:, 1:n])
    σ::Vector{Int}                   # photon order permutation, length n
    out::Vector{Int}                 # sampled output modes, length n
    perm_full::Matrix{ComplexF64}    # (n−1) × n scratch holding A[out[1:k−1], σ[1:k]]
    sub_buf::Matrix{ComplexF64}      # (n−1) × (n−1) scratch for ryser input
    v_perms::Vector{ComplexF64}      # length n: per of (k−1)×(k−1) submatrices
    weights::Vector{Float64}         # length m: per-mode unnormalized weight
end

function CCSamplerWorkspace(U_in::AbstractMatrix, n::Int)
    m = size(U_in, 1)
    @assert size(U_in, 2) >= n
    A = ComplexF64.(U_in[:, 1:n])
    return CCSamplerWorkspace(A,
        Vector{Int}(undef, n), Vector{Int}(undef, n),
        Matrix{ComplexF64}(undef, max(n - 1, 1), n),
        Matrix{ComplexF64}(undef, max(n - 1, 1), max(n - 1, 1)),
        Vector{ComplexF64}(undef, n),
        Vector{Float64}(undef, m))
end

# Inverse-CDF sample of an integer in 1..m proportional to the (positive) weights.
function _wsample_cdf(weights::AbstractVector{Float64}, m::Int)
    total = 0.0
    @inbounds @simd for i in 1:m; total += weights[i]; end
    u = rand() * total
    cum = 0.0
    @inbounds for i in 1:m
        cum += weights[i]
        u <= cum && return i
    end
    return m
end

# Sample n output modes in-place into ws.out (returns sorted ws.out).
function cc_sample!(ws::CCSamplerWorkspace)
    m, n = size(ws.A)
    σ = ws.σ
    @inbounds for i in 1:n; σ[i] = i; end
    Random.shuffle!(σ)

    # First photon (input mode σ[1]): weight ∝ |U[i, σ[1]]|².
    s1 = σ[1]
    @inbounds @simd for i in 1:m
        ws.weights[i] = abs2(ws.A[i, s1])
    end
    ws.out[1] = _wsample_cdf(ws.weights, m)

    # Subsequent photons via the chain-rule expansion of |Per|².
    @inbounds for k in 2:n
        # Build (k−1)×k matrix: rows = previously placed photons, cols = σ[1..k]
        for i in 1:(k - 1)
            row = ws.out[i]
            for j in 1:k
                ws.perm_full[i, j] = ws.A[row, σ[j]]
            end
        end
        # k permanents of (k−1)×(k−1) submatrices (one column removed).
        if k == 2
            # (k−1)×(k−1) is 1×1 → permanent is just the entry.
            ws.v_perms[1] = ws.perm_full[1, 2]
            ws.v_perms[2] = ws.perm_full[1, 1]
        else
            for skip in 1:k
                jj = 1
                for j in 1:k
                    j == skip && continue
                    for i in 1:(k - 1)
                        ws.sub_buf[i, jj] = ws.perm_full[i, j]
                    end
                    jj += 1
                end
                ws.v_perms[skip] = ryser(@view ws.sub_buf[1:(k - 1), 1:(k - 1)])
            end
        end
        # weights[i] = |Σ_j U[i, σ[j]] · v_perms[j]|²
        for i in 1:m
            s = ComplexF64(0)
            for j in 1:k
                s += ws.A[i, σ[j]] * ws.v_perms[j]
            end
            ws.weights[i] = abs2(s)
        end
        ws.out[k] = _wsample_cdf(ws.weights, m)
    end

    sort!(ws.out)
    return ws.out
end

# Legacy interface (for testing).
function estimate_S(U::AbstractMatrix, omega::AbstractVector, N_fourier::Integer,
                    x0::Integer, n_photons::Int, M::Int)
    total = complex(0.0)
    for _ in 1:M
        total += Z_sample(U, omega, N_fourier, x0, n_photons)
    end
    return real(total / M)
end

# ─── Binning (Section 1.1) ───────────────────────────────────────────
# Divide total_states into d equal-width bins.
function bin_edges(total_states::Integer, d::Int)
    T = typeof(total_states)
    base_width = total_states ÷ d
    remainder = total_states % d

    edges = Vector{T}(undef, d + 1)
    edges[1] = zero(T)
    for j in 1:d
        w = base_width + (j <= remainder ? oneunit(T) : zero(T))
        edges[j + 1] = edges[j] + w
    end
    return edges
end

# ─── Unranking: find the i-th state in f-value order ─────────────────
# For base-(n+1) encoding: f(s) = sum_j s_j * (n+1)^(j-1).
# Returns the occupation vector s of the rank-th state (0-indexed)
# among all weak compositions of n_photons into m non-negative parts,
# ordered by f-value (most significant digit = last mode).
function unrank_composition(rank::Integer, n_photons::Int, m::Int)
    s = zeros(Int, m)
    remaining_n = n_photons
    remaining_rank = big(rank)
    for pos in m:-1:2
        for k in 0:remaining_n
            count = binomial(big(remaining_n - k + pos - 2), big(pos - 2))
            if remaining_rank < count
                s[pos] = k
                remaining_n -= k
                break
            end
            remaining_rank -= count
        end
    end
    s[1] = remaining_n
    return s
end

# f-value of state s under base encoding: f(s) = sum_j s_j * base^(j-1).
# Returns the same integer type as base (Int64, Int128, or BigInt).
function f_value(s::AbstractVector, base::Integer)
    T = typeof(base)
    result = zero(T)
    pow = oneunit(T)
    for j in eachindex(s)
        result += T(s[j]) * pow
        pow *= base
    end
    return result
end

# ─── Auto-select integer type for N ──────────────────────────────────
# Picks the smallest integer type that fits n * (n+1)^(m-1) + 1.
function _compute_N(n::Int, base::Int, m::Int)
    # F64-E: Float64 only for magnitude comparison; thresholds have >2× safety margin.
    Nf = Float64(n) * Float64(base)^(m - 1)
    if Nf < 4.0e18
        return n * base^(m - 1) + 1
    elseif Nf < 1.7e38
        return Int128(n) * Int128(base)^(m - 1) + Int128(1)
    else
        return BigInt(n) * BigInt(base)^(m - 1) + BigInt(1)
    end
end

# ─── Full algorithm: find the most probable bin (Corollary 1) ─────────
function find_most_probable_bin(U::AbstractMatrix, m::Int, n::Int, d::Int;
                                M::Int=10000)
    base = n + 1
    # _compute_N returns Int64 / Int128 / BigInt — use a function barrier
    # so that Julia fully specializes the inner loop for the concrete type.
    N_fourier = _compute_N(n, base, m)
    return _find_most_probable_bin(U, m, n, d, M, N_fourier)
end

# Function barrier: fully specialized for concrete integer type T.
function _find_most_probable_bin(U::AbstractMatrix, m::Int, n::Int, d::Int, M::Int,
                                 N_fourier::T) where T <: Integer
    base_typed = T(n + 1)

    # Total number of output states (stars-and-bars): C(n+m-1, n)
    total_states = binomial(big(n + m - 1), big(n))
    edges = bin_edges(total_states, d)

    ctx = SamplingContext(N_fourier)
    U_in = U[:, 1:n]

    # Estimate CDF S(x_j) at each bin boundary via unranking
    S_values = Vector{Float64}(undef, d + 1)
    S_values[1] = 0.0

    for j in 1:d
        if edges[j + 1] == 0
            S_values[j + 1] = 0.0
        elseif edges[j + 1] == total_states
            S_values[j + 1] = 1.0
        else
            s = unrank_composition(edges[j + 1] - 1, n, m)
            x_j = f_value(s, base_typed)
            S_values[j + 1] = estimate_S(U_in, base_typed, N_fourier, x_j, n, M, ctx)
        end
    end

    bin_probs = [S_values[j + 1] - S_values[j] for j in 1:d]
    j_star = argmax(bin_probs)
    return j_star, bin_probs
end

# ─── State enumeration (for testing / small instances) ────────────────

# Enumerate all non-collision states (at most 1 photon per mode).
function _enumerate_noncollision_states(m::Int, n::Int)
    states = Vector{Vector{Int}}()
    _enumerate_recursive!(states, Int[], m, n, 1)
    return states
end

function _enumerate_recursive!(states, current, m, n, start)
    if length(current) == m
        if sum(current) == n
            push!(states, copy(current))
        end
        return
    end
    remaining_modes = m - length(current)
    remaining_photons = n - sum(current; init=0)
    # Pruning
    if remaining_photons < 0 || remaining_photons > remaining_modes
        return
    end
    for s in 0:min(1, remaining_photons)
        push!(current, s)
        _enumerate_recursive!(states, current, m, n, start)
        pop!(current)
    end
end

# ─── Time estimation ─────────────────────────────────────────────────

# Measure per-sample cost (seconds) via a short pilot run.
function _pilot_time_per_sample(m::Int, n::Int; M_pilot::Int=5000)
    base = n + 1
    N = _compute_N(n, base, m)
    return _pilot_inner(m, n, N, M_pilot)
end

function _pilot_inner(m::Int, n::Int, N::T, M_pilot::Int) where T <: Integer
    base_typed = T(n + 1)
    ctx = SamplingContext(N)
    U_in = RandHaar(m).U[:, 1:n]
    return @elapsed(estimate_S(U_in, base_typed, N, N ÷ 2, n, M_pilot, ctx)) / M_pilot
end

# Estimate wall time (seconds) for find_most_probable_bin with a fixed M.
function estimate_time(m::Int, n::Int, d::Int, M::Int; M_pilot::Int=5000)
    t_per_sample = _pilot_time_per_sample(m, n; M_pilot)
    return (d - 1) * M * t_per_sample
end

# Estimate wall time (seconds) to reach a target max bin-probability error ε.
# Runs n_trials pilot pipeline calls at M_pilot to calibrate error, then
# extrapolates M_required via the 1/√M scaling law.
function estimate_time_to_precision(m::Int, n::Int, d::Int, ε::Float64;
                                    M_pilot::Int=5000, n_trials::Int=5)
    t_per_sample = _pilot_time_per_sample(m, n; M_pilot)

    # Run pilot pipelines to estimate typical error at M_pilot
    U = RandHaar(m).U
    pilot_errs = Float64[]
    for _ in 1:n_trials
        _, bp = find_most_probable_bin(U, m, n, d; M=M_pilot)
        push!(pilot_errs, std(bp))
    end
    err_pilot = mean(pilot_errs)

    M_required = ceil(Int, M_pilot * (err_pilot / ε)^2)
    M_required = max(M_required, M_pilot)
    t_total = (d - 1) * M_required * t_per_sample
    return (; time=t_total, M=M_required)
end
