using BosonSampling
using BenchmarkTools
using BenchmarkPlots, StatsPlots
using JLD

push!(LOAD_PATH, "./benchmarks")
run(`python ./benchmarks/permanents_benchmarks/thewalrus_perm.py`)
run(`python ./benchmarks/permanents_benchmarks/pcvl_perm.py`)
run(`python ./benchmarks/permanents_benchmarks/compute_perm.py`)

BenchmarkTools.DEFAULT_PARAMETERS.samples = 1000
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1
BenchmarkTools.DEFAULT_PARAMETERS.evals = 5

suite_permanent = BenchmarkGroup(["string"])
suite_permanent["ryser-jl"] = BenchmarkGroup(["string"])
suite_permanent["glynn-jl"] = BenchmarkGroup(["string"])
for n in 1:30
    A = RandHaar(n)
    suite_permanent["ryser-jl"]["n=$n"] = @benchmarkable ryser($A.U)
    suite_permanent["glynn-jl"]["n=$n"] = @benchmarkable fast_glynn_perm($A.U)
end
res = run(suite_permanent, verbose=true, seconds=1)
save("benchmarks/data/permanents_bench.jld", "permanents_benchmarkGroup",res)
