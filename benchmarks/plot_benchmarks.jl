using LaTeXStrings
using JLD
using BenchmarkTools
using Plots
push!(LOAD_PATH, "./benchmarks")

d = load("benchmarks/data/permanents_bench.jld")
ryser_jl = [mean(d["permanents_benchmarkGroup"]["ryser-jl"]["n=$n"]).time * 10^(-9) for n in 1:30]
glynn_jl = [mean(d["permanents_benchmarkGroup"]["glynn-jl"]["n=$n"]).time * 10^(-9) for n in 1:30]

f = open("benchmarks/data/data_python.txt")
f = readlines(f)
python_data = map(e -> parse(Float64,e), f)
python_data = vcat(zeros(14), python_data)

f = open("benchmarks/data/data_matlab.txt")
f = readlines(f)
matlab_data = map(e -> parse(Float64,e), f)
matlab_data = vcat(zeros(14), matlab_data)

fig = plot(xlabel=L"Matrix size $n$",
    ylabel=L"$log_{10}($Time$)$ (s)",
    legend=:topleft,
    yaxis =:log10,
    xlims = (14.5,25.5),
    xticks = 14:1:26,
    ylim=[10^(-4),600],
    # yticks = 10^(-4):600,
    yminorticks = 5,
    yminorgrid=true,
    dpi=800)

plot!(ryser_jl, markershape=:+, label="ryser-jl", dpi=800)
scatter!(python_data, label="ryser-py", dpi=800)
scatter!(matlab_data, label="ryser-m", dpi=800)
savefig("docs/publication/package/images/permanents_bench1.png")

f = open("benchmarks/data/thewalrus_data.txt")
f = readlines(f)
the_walrus_data = map(e->parse(Float64,e), f)

f = open("benchmarks/data/glynn-pcvl.txt")
f = readlines(f)
glynn_pcvl = map(e->parse(Float64,e), f)

f = open("benchmarks/data/ryser4-pcvl.txt")
f = readlines(f)
ryser4_pcvl = map(e->parse(Float64,e), f)

fig = plot(ylabel=L"$log_{10}($Time$)$ (s)",
    legend=:topleft,
    yaxis =:log10,
    xlims = (0.5,30.5),
    xticks = 0:5:30,
    ylim=[10^(-7),600],
    # yticks = 10^(-4):600,
    yminorticks = 5,
    yminorgrid=true,
    dpi=800)

plot!(ryser_jl, markershape=:+, label="ryser-jl", dpi=800)
scatter!(the_walrus_data, markershape=:utriangle, marker=3, label="thewalrus", markerstrokewidth=0, dpi=800)
scatter!(glynn_pcvl,markershape=:star4, marker=3, label="glynn-pcvl", markerstrokewidth=0, dpi=800)
scatter!(ryser4_pcvl, marker=3, label="ryser-4-pcvl", markerstrokewidth=0, dpi=800)

d = load("benchmarks/data/samplers_1.jld")
data_1 = [mean(d["sampler1"]["n=$n"]).time * 10^(-9) for n in 1:30]
d = load("benchmarks/data/samplers_4.jld")
data_4 = [mean(d["sampler4"]["n=$n"]).time * 10^(-9) for n in 1:30]
d = load("benchmarks/data/samplers_8.jld")
data_8 = [mean(d["sampler8"]["n=$n"]).time * 10^(-9) for n in 1:30]

f = open("benchmarks/data/cliffords-pcvl.txt")
f = readlines(f)
cliffords_pcvl = map(e->parse(Float64,e), f)

fig_samp = plot(xlabel=L"Matrix size $n$",
    ylabel=L"$log_{10}($Time$)$ (s)",
    legend=:topleft,
    yaxis =:log10,
    xlims = (0.5,30.5),
    xticks = 0:5:30,
    ylim=[10^(-7),600],
    # yticks = 10^(-4):600,
    yminorticks = 5,
    yminorgrid=true,
    dpi=800)

plot!(data_1, markershape=:+, label="clifford's-jl-1", dpi=800)
plot!(data_4, markershape=:+, label="clifford's-jl-4", dpi=800)
plot!(data_8, markershape=:+, label="clifford's-jl-8", dpi=800)
scatter!(cliffords_pcvl, label="clifford's-pcvl", markerstrokewidth=0, markersize=3, dpi=800)

plot(fig, fig_samp, layout= grid(2,1, heights=[0.5,0.5], widths=[0.5,0.5]), title=["(a)" "(b)"])
plot!(size=(850,900))
savefig("docs/publication/package/images/exist_soft.png")
