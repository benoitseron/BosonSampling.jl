begin
    using Revise

    using BosonSampling
    using Plots
    using ProgressMeter
    using Distributions
    using Random
    using Test
    using ArgCheck
    using StatsBase
    using ColorSchemes
    using Interpolations
    using Dierckx
    using LinearAlgebra
    using PrettyTables
    using LaTeXStrings
    using JLD
    using AutoHashEquals
    using LinearRegression
    using DataStructures
    using Parameters
    using UnPack

end

n_array = collect(3:16)
for input_type in [Distinguishable, Bosonic]

    function get_sample_time(n, niter = 100)

        m = n

        (@elapsed for i in 1:niter get_sample_loop(LoopSamplingParameters(n=n, input_type = input_type, η_loss_bs = 0.9 .* ones(m-1), η_loss_lines = 0.9 .* ones(m))) end)/niter

    end

    get_sample_time(3) # compilation


    time_array = [get_sample_time(n) for n in n_array]

    plot(n_array, time_array, yaxis = :log10, label = false)
    xaxis!("n_photons (with loss)")
    yaxis!("time one sample (s)")
    title!("$(input_type)")
    savefig("benchmarks/sampler_benchmarks/loop/images/sample_time_lossy_$(input_type).png")

end
