# Publication-clean variance scaling: Var(Z) vs n (m=n), averaged over several
# Haar instances per n, fit to a power law, pushed to n=200.  Confirms that the
# samples-per-fixed-error  M = Var(Z)/ε²  grows POLYNOMIALLY, so one value of the
# cumulative function S(x0) at fixed error is polynomial-time at all sizes.
#
# Run from the BosonSampling.jl package root:
#   JULIA_NUM_THREADS=N julia --project=. docs/publication/one_way_function/variance_scaling_clean.jl

using LinearAlgebra, Random, Printf, Statistics
using BosonSampling
using Plots

cd(@__DIR__)   # write figures next to this script, regardless of launch dir

default(size=(720,460), bottom_margin=6Plots.mm, left_margin=9Plots.mm,
        right_margin=4Plots.mm, top_margin=4Plots.mm,
        titlefont=font("sans-serif",11), guidefont=font("sans-serif",10),
        tickfont=font("sans-serif",9), legendfont=font("sans-serif",8),
        framestyle=:box, grid=true, gridalpha=0.25)

function setup(n, m, seed)
    Random.seed!(seed)
    base = n + 1; N = _compute_N(n, base, m); T = typeof(N)
    U_in = Matrix{ComplexF64}(RandHaar(m).U[:, 1:n])
    return (U_in, T(base), N, n, SamplingContext(N))
end

# x0 fraction where S≈0.5 (worst MC case), coarse log scan.
function find_median_q(U_in, base, N, nn, ctx; M=30_000)
    T = typeof(N); bestq, bestd = 0.5, Inf
    for q in 10.0 .^ range(-7, log10(0.5), length=16)
        x0 = T(round(BigInt, BigInt(N) * q))
        S = estimate_S(U_in, base, N, x0, nn, M, ctx)
        abs(S-0.5) < bestd && ((bestd, bestq) = (abs(S-0.5), q))
    end
    bestq
end

function var_at(U_in, base, N, nn, ctx, q; M_var=150_000, trials=20)
    T = typeof(N); x0 = T(round(BigInt, BigInt(N)*q))
    estimate_S(U_in, base, N, x0, nn, 2000, ctx)
    Sv = [estimate_S(U_in, base, N, x0, nn, M_var, ctx) for _ in 1:trials]
    var(Sv)*M_var
end

function main()
    println("threads = $(Threads.nthreads())\n")
    # span the Float64 boundary (N crosses 1.8e308 at n≈143): 120,140 use the finite
    # path; 150,160 use the Inf path — continuity here validates the huge-N branches.
    ns = [20,40,60,80,100,120,140,150,160,200,250,300]
    n_inst = 4
    meanvar = Float64[]; stdvar = Float64[]
    println("    n   N-type    <Var(Z)>   (per-instance)            M_req(1%)")
    for n in ns
        vs = Float64[]
        Nt = ""
        for s in 1:n_inst
            U_in, base, N, nn, ctx = setup(n, n, 1000 + 17*n + s)
            Nt = string(typeof(N))
            q = find_median_q(U_in, base, N, nn, ctx)
            push!(vs, var_at(U_in, base, N, nn, ctx, q))
        end
        mv = mean(vs); push!(meanvar, mv); push!(stdvar, std(vs))
        @printf("  %4d  %8s   %8.1f   [%s]   %.2e\n",
                n, Nt, mv, join((@sprintf("%.0f",x) for x in vs), " "), mv/1e-4)
    end

    # power-law fit  Var ≈ a·n^p  (log-log least squares)
    ln = log.(ns); lv = log.(meanvar); k = length(ns)
    p = (k*sum(ln.*lv)-sum(ln)*sum(lv))/(k*sum(ln.^2)-sum(ln)^2)
    a = exp(sum(lv)/k - p*sum(ln)/k)
    ss_res = sum((lv .- (log(a) .+ p.*ln)).^2)
    ss_tot = sum((lv .- mean(lv)).^2)
    R2 = 1 - ss_res/ss_tot
    @printf("\nFIT:  Var(Z) ≈ %.4f · n^%.3f      (R² = %.4f)\n", a, p, R2)
    @printf("  => M_required(ε) ≈ %.3f · n^%.2f / ε²   samples\n", a, p)

    # plot: measured points + fit, log-log, with exp reference
    nf = range(ns[1], ns[end], length=200)
    p1 = scatter(ns, meanvar, yerror=stdvar, ms=5, color=:steelblue,
                 label="measured ⟨Var(Z)⟩  (4 Haar instances)", legend=:topleft,
                 xscale=:log10, yscale=:log10, xlabel="n  (photons = modes)",
                 ylabel="Var(Z)  at  S≈0.5",
                 title="Estimator variance is polynomial in n")
    plot!(p1, nf, a .* nf.^p, lw=2, color=:tomato, ls=:dash,
          label=@sprintf("fit  %.3f·n^{%.2f}  (R²=%.3f)", a, p, R2))
    plot!(p1, nf, meanvar[1] .* 2 .^((nf .- ns[1])./15), lw=1.5, color=:gray, ls=:dot,
          label="exponential ref  2^{n/15}")
    savefig(p1, "variance_scaling_clean.png")
    savefig(p1, "variance_scaling_clean.pdf")
    println("\nSaved variance_scaling_clean.png / .pdf")

    # headline at n=m=100 and 200 (per-sample time measured)
    for n in (100, 200, 300)
        U_in, base, N, nn, ctx = setup(n, n, 999)
        q = find_median_q(U_in, base, N, nn, ctx)
        x0 = typeof(N)(round(BigInt, BigInt(N)*q))
        estimate_S(U_in, base, N, x0, nn, 2000, ctx)
        t = @elapsed estimate_S(U_in, base, N, x0, nn, 50_000, ctx); tps = t/50_000
        Vz = a*n^p; Mreq = Vz/1e-4; tval = Mreq*tps
        @printf("  n=m=%d:  Var≈%.0f  M_req(1%%)≈%.2e  t/samp=%.1fµs  =>  one S(x0) ≈ %s on %d threads\n",
                n, Vz, Mreq, tps*1e6,
                tval<60 ? @sprintf("%.1f s",tval) : @sprintf("%.1f min",tval/60),
                Threads.nthreads())
    end
end

main()
