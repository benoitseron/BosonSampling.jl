using BosonSampling
using BenchmarkTools
using BenchmarkPlots, StatsPlots
using Permanents
using LaTeXStrings
using JLD

push!(LOAD_PATH, "./benchmarks")
run(`python ./benchmarks/thewalrus_data.py`)
run(`python ./benchmarks/pcvl_data.py`)

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
save("benchmarks/permanents_bench.jld", "permanents_benchmarkGroup",res)

d = load("benchmarks/permanents_bench.jld")
ryser_jl = [mean(d["permanents_benchmarkGroup"]["ryser-jl"]["n=$n"]).time * 10^(-9) for n in 1:30]
glynn_jl = [mean(d["permanents_benchmarkGroup"]["glynn-jl"]["n=$n"]).time * 10^(-9) for n in 1:30]

f = open("benchmarks/thewalrus_data.txt")
f = readlines(f)
the_walrus_data = map(e->parse(Float64,e), f)

f = open("benchmarks/glynn-pcvl.txt")
f = readlines(f)
glynn_pcvl = map(e->parse(Float64,e), f)

f = open("benchmarks/ryser4-pcvl.txt")
f = readlines(f)
ryser4_pcvl = map(e->parse(Float64,e), f)


fig = plot(xlabel="Matrix size n", ylabel="log(Time) (s)", yaxis=:log, legend=:bottomright, dpi=300)
scatter!(the_walrus_data, marker=3, label="thewalrus", dpi=300)
scatter!(ryser_jl, marker=3, label="ryser-jl", dpi=300)
scatter!(glynn_jl, marker=3, label="glynn-jl", dpi=300)
scatter!(glynn_pcvl, marker=3, label="glynn-pcvl", dpi=300)
scatter!(ryser4_pcvl, marker=3, label="ryser-4-pcvl", dpi=300)
savefig(fig, "docs/src/benchmarks/bench_perm.png")

# f = open("data_python.txt")
# f = readlines(f)
# python_data = split(f[1], " ")
# python_data = map(e->parse(Float64,e), python_data)
#
# f = open("data_matlab.txt")
# f = readlines(f)
# matlab_data = split(f[1], " ")
# matlab_data = map(e->parse(Float64,e), matlab_data)
#
# perm_fig = plot(xlabel="n", ylabel="time (s)", legend=:topleft, dpi=300)
# scatter!(15:25, data_ryser, label="Julia", dpi=300)
# scatter!(15:25, python_data, label="python", dpi=300)
# scatter!(15:25, matlab_data, label="Matlab", dpi=300)
# savefig(perm_fig, "docs/src/benchmarks/compute_perm.png")

suite_sampler = BenchmarkGroup()
suite_sampler["cliffords_sampler"] = BenchmarkGroup(["string"])
suite_sampler["noisy_sampler"] = BenchmarkGroup(["string"])

for n in 1:30
    interf = Fourier(n)
    input = Input{Bosonic}(first_modes(n,n))
    suite_sampler["cliffords_sampler"]["n=$n"] = @benchmarkable cliffords_sampler(input=$input, interf=$interf)
    input = Input{OneParameterInterpolation}(first_modes(n,n), 0.976)
    suite_sampler["noisy_sampler"]["n=$n,reflec=0.755"] = @benchmarkable noisy_sampler(input=$input, reflectivity=0.755, interf=$interf)
end

res = run(suite_sampler, verbose=true, seconds = 1)
data_bosonic = []
data_noisy = []

for j in 1:30
    mean_bosonic = mean(res["cliffords_sampler"]["n=$j"])
    push!(data_bosonic, (mean_bosonic.time)/10^9)
    mean_noisy = mean(res["noisy_sampler"]["n=$j,reflec=0.755"])
    push!(data_noisy, (mean_noisy.time)/10^9)
end

save("benchmarks/samplers.jld", "samplers", res)

d = load("benchmarks/samplers.jld")
data_bosonic = [mean(d["samplers"]["cliffords_sampler"]["n=$n"]).time * 10^(-9) for n in 1:30]
data_noisy = [mean(d["samplers"]["noisy_sampler"]["n=$n,reflec=0.755"]).time * 10^(-9) for n in 1:30]

fig_samp = plot(xlabel="n", ylabel="time (s)", yaxis=:log, legend=:topleft, dpi=300)
scatter!(2:30, data_bosonic, label="cliffords sampler", dpi=300, marker=2)
scatter!(2:30, data_noisy, label=L"noisy sampler, $η=0.755", dpi=300, marker=2)
savefig(fig_samp, "docs/src/benchmarks/sampler.png")
