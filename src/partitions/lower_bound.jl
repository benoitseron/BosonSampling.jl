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


n = 10
m = 12
x = 0.9


function get_avg_proba(n,m,x, niter = 1000)

    results = zeros(niter)

    for j in 1:niter
        i = Input{OneParameterInterpolation}(first_modes(n,m),x)

        s1 = Subset(first_modes(Int(m/2), m))
        interf = RandHaar(m)
        part = Partition(s1)
        o = PartitionCount(PartitionOccupancy(ModeOccupation([Int(n/2)]),  n, part))
        ev = Event(i,o,interf)

        compute_probability!(ev)

        results[j] = ev.proba_params.probability
    end

    mean(results)

end

x_array = collect(range(0.9,1, length = 100))
proba_array = get_avg_proba.(n,m,x_array)

plot(x_array, proba_array)

n = 10

m_array = collect(n:4:10n)

# x = 0.8
# proba_array = get_avg_proba.(n,m_array,x)
# plot(m_array, proba_array, label = "x = $x")
# x = 1
# proba_array = get_avg_proba.(n,m_array,x)
# plot!(m_array, proba_array, label = "x = $x")

x = 0.9
tvd_lower_bound_array =  abs.(get_avg_proba.(n,m_array,x) - get_avg_proba.(n,m_array,1))

plot(m_array, tvd_lower_bound_array)

n_array = collect(2:2:16)

plot(n_array, get_avg_proba.(n_array,n_array,1))
