using BosonSampling
using BenchmarkTools
using BenchmarkPlots, StatsPlots
using JLD

push!(LOAD_PATH, "./benchmarks/")
run(`python ./benchmarks/sampler_benchmarks/sampler_pcvl.py`)

BenchmarkTools.DEFAULT_PARAMETERS.samples = 1000
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1
BenchmarkTools.DEFAULT_PARAMETERS.evals = 5

suite_sampler1 = BenchmarkGroup()
Threads.nthreads() = 1
for n in 1:30
    interf = RandHaar(n)
    input = Input{Bosonic}(first_modes(n,n))
    suite_sampler1["n=$n"] = @benchmarkable cliffords_sampler(input=$input, interf=$interf)
end
res = run(suite_sampler1, verbose=true, seconds = 1)
save("benchmarks/data/samplers_1.jld", "sampler1",res)

suite_sampler4 = BenchmarkGroup()
Threads.nthreads() = 4
for n in 1:30
    interf = RandHaar(n)
    input = Input{Bosonic}(first_modes(n,n))
    suite_sampler4["n=$n"] = @benchmarkable cliffords_sampler(input=$input, interf=$interf)
end
res = run(suite_sampler4, verbose=true, seconds = 1)
save("benchmarks/data/samplers_4.jld", "sampler4",res)

suite_sampler8 = BenchmarkGroup()
Threads.nthreads() = 8
for n in 1:30
    interf = RandHaar(n)
    input = Input{Bosonic}(first_modes(n,n))
    suite_sampler8["n=$n"] = @benchmarkable cliffords_sampler(input=$input, interf=$interf)
end
res = run(suite_sampler8, verbose=true, seconds = 1)
save("benchmarks/data/samplers_8.jld", "sampler8", res)
