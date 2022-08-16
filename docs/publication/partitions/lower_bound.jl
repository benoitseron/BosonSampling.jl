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
m = 20

###### checking if the 1 distinguishable formula holds ######

# set the interf before the rest
interf = RandHaar(m)

# S matrix with one dist photon
"""
    gram_one_dist_photon(dist_photon, n)

Gram matrix with all indistinguishable photons except the one at position `dist_photon`.
"""
function gram_one_dist_photon(dist_photon, n)

    mat = ones(ComplexF64,(n,n))
    mat[dist_photon, :] .= 0
    mat[:, dist_photon] .= 0
    mat[dist_photon, dist_photon] = 1

    mat

end

ib = Input{Bosonic}(first_modes(n,m))
i(dist_photon) = Input{UserDefinedGramMatrix}(first_modes(n,m), gram_one_dist_photon(dist_photon, n))

s1 = Subset(first_modes(Int(m/2), m))

part = Partition(s1)
o = PartitionCount(PartitionOccupancy(ModeOccupation([Int(n/2)]),  n, part))


"""
Exact probability with `x`.
"""
function actual_probability(x)

    ix = Input{OneParameterInterpolation}(first_modes(n,m),x)
    ev = Event(ix,o,interf)

    compute_probability!(ev)
    ev.proba_params.probability

end

"""
Approximate probability proposal.
"""
function approximate_probability(x)

    ϵ = 1 - x

    ib = Input{Bosonic}(first_modes(n,m))

    i(dist_photon) = Input{UserDefinedGramMatrix}(first_modes(n,m), gram_one_dist_photon(dist_photon, n))

    evb = Event(ib,o,interf)
    compute_probability!(evb)
    pb = evb.proba_params.probability

    p_dist_array = zeros(n)
    for dist_photon in 1:n

        evd = Event(i(dist_photon),o,interf)
        compute_probability!(evd)
        p_dist_array[dist_photon] = evd.proba_params.probability

    end

    p_dist = mean(p_dist_array)

    pb * (1 - n*ϵ) + n*ϵ * p_dist

end

x = 0.9
actual_probability(x)
approximate_probability(x)

x_array = collect(range(0.5, 1, length = 30))
plot(x_array, actual_probability.(x_array), label = "p(x)")
plot!(x_array, approximate_probability.(x_array), label = "papprox(x)")
title!("validity of Eq. 157, n = $n, m = $m")
xlabel!("x")
savefig("src/partitions/images/validity_approximate_partition_formula_n = $n, m = $m.png")



###### lower bound ######

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
