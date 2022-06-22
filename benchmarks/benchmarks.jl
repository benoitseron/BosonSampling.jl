using BosonSampling
using BenchmarkTools
using BenchmarkPlots, StatsPlots
using Permanents
using LaTeXStrings
using JLD
using FileIO

BenchmarkTools.DEFAULT_PARAMETERS.samples = 50
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1
BenchmarkTools.DEFAULT_PARAMETERS.evals = 1

suite_permanent = BenchmarkGroup(["string"])
for n in 15:25
    A = RandHaar(n)
    suite_permanent["n=$n"] = @benchmarkable ryser($A.U)
end

res = run(suite_permanent, verbose=true, seconds=1)
data_mean = []

for i in 15:25
    mean_res = mean(res["n=$i"])
    push!(data_mean, (mean_res.time)/10^9)
end

f = open("data_python.txt")
f = readlines(f)
python_data = split(f[1], " ")
python_data = map(e->parse(Float64,e), python_data)

f = open("data_matlab.txt")
f = readlines(f)
matlab_data = split(f[1], " ")
matlab_data = map(e->parse(Float64,e), matlab_data)

perm_fig = plot(xlabel="n", ylabel="time (s)", legend=:topleft, dpi=300)
scatter!(15:25, data_mean, label="Julia", dpi=300)
scatter!(15:25, python_data, label="python", dpi=300)
scatter!(15:25, matlab_data, label="Matlab", dpi=300)
savefig(perm_fig, "docs/src/benchmarks/compute_perm.png")


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

jldopen("benchmarks/data_bosonic.jld", "w") do file
    write(file, "data_bosonic", data_bosonic)
end

d = jldopen("benchmarks/data_bosonic.jld", "r") do file
    read(file, "data_bosonic")
end

fig_samp = plot(xlabel="n", ylabel="time (s)", legend=:topleft, dpi=300)
scatter!(2:30, data_bosonic, label="cliffords sampler", dpi=300)
scatter!(2:30, data_noisy, label=L"noisy sampler, $η=0.755", dpi=300)
savefig(fig_samp, "docs/src/benchmarks/sampler.png")
