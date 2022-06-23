using BosonSampling
using BenchmarkTools
using BenchmarkPlots, StatsPlots
using Permanents
using LaTeXStrings
using JLD

BenchmarkTools.DEFAULT_PARAMETERS.samples = 1000
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5
BenchmarkTools.DEFAULT_PARAMETERS.evals = 10

suite_permanent = BenchmarkGroup(["string"])
suite_permanent["ryser-jl"] = BenchmarkGroup(["string"])
suite_permanent["glynn-jl"] = BenchmarkGroup(["string"])
for n in 2:30
    A = RandHaar(n)
    suite_permanent["ryser-jl"]["n=$n"] = @benchmarkable ryser($A.U)
    suite_permanent["glynn-jl"]["n=$n"] = @benchmarkable fast_glynn_perm($A.U)
end

res = run(suite_permanent, verbose=true, seconds=1)
data_ryser = []
data_glynn = []

# for i in 2:30
#     mean_res = mean(res["ryser-jl"]["n=$i"])
#     push!(data_ryser, (mean_res.time)/10^9)
#     mean_res = mean(res["glynn-jl"]["n=$i"])
#     push!(data_glynn, (mean_res.time)/10^9)
# end

save("benchmarks/permanents_bench.jld", "permanents_benchmarkGroup",res)
d = load("benchmarks/permanents_bench.jld")
data_ryser = [mean(d["permanents_benchmarkGroup"]["ryser-jl"]["n=$n"]).time * 10^(-9) for n in 2:30]

f = open("benchmarks/data_thewalrus.txt")
f = readlines(f)
the_walrus_data = split(f[:], " ")
the_walrus_data = map(e->parse(Float64,e), the_walrus_data)

f = open("benchmarks/data_pcvl.txt")
f = readlines(f)
pcvl = split(f[:], " ")
pcvl = map(e->parse(Float64,e), pcvl)

data_ryser = map(e -> log(e), data_ryser)
the_walrus_data = map(e -> log(e), the_walrus_data)
pcvl = map(e -> log(e), pcvl)

fig = plot(xlabel="Matrix size n", ylabel="log(Time) (s)", legend=:bottomright)
scatter!(the_walrus_data, marker=2)
scatter!(data_ryser, marker=2)
scatter!(pcvl, marker=2)

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

for n in 2:30
    interf = Fourier(n)
    input = Input{Bosonic}(first_modes(n,n))
    suite_sampler["cliffords_sampler"]["n=$n"] = @benchmarkable cliffords_sampler(input=$input, interf=$interf)
    input = Input{OneParameterInterpolation}(first_modes(n,n), 0.976)
    suite_sampler["noisy_sampler"]["n=$n,reflec=0.755"] = @benchmarkable noisy_sampler(input=$input, reflectivity=0.755, interf=$interf)
end

res = run(suite_sampler, verbose=true, seconds = 1)
data_bosonic = []
data_noisy = []

for j in 2:30
    mean_bosonic = mean(res["cliffords_sampler"]["n=$j"])
    push!(data_bosonic, (mean_bosonic.time)/10^9)
    mean_noisy = mean(res["noisy_sampler"]["n=$j,reflec=0.755"])
    push!(data_noisy, (mean_noisy.time)/10^9)
end

save("benchmarks/samplers.jld", "samplers", res)

fig_samp = plot(xlabel="n", ylabel="time (s)", legend=:topleft, dpi=300)
scatter!(2:30, data_bosonic, label="cliffords sampler", dpi=300)
scatter!(2:30, data_noisy, label=L"noisy sampler, $η=0.755", dpi=300)
savefig(fig_samp, "docs/src/benchmarks/sampler.png")
