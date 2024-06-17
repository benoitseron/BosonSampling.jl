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





function tvd_two_partitions(n,m,n_subsets)
    ib = Input{Bosonic}(first_modes(n,m))
    interf = RandHaar(m)
    part_a = random_partition(m,n_subsets)
    part_b = part_a
    while part_a == part_b
        part_b = random_partition(m,n_subsets)
    end

    o_a = PartitionCountsAll(part_a)
    o_b = PartitionCountsAll(part_b)

    ev_a = Event(ib,o_a,interf)
    ev_b = Event(ib,o_b,interf)

    compute_probability!(ev_a)
    compute_probability!(ev_b)

    tvd(ev_a.proba_params.probability.proba, ev_b.proba_params.probability.proba)

end

function variation_two_partitions(n,m, n_subsets, n_iter)
    tvd_array = zeros(n_iter)
    for i in 1:n_iter
        tvd_array[i] = tvd_two_partitions(n,m,n_subsets)
    end

    mean_tvd = mean(tvd_array)
    std_tvd = std(tvd_array)
    return mean_tvd, std_tvd
end

n_array = 2:1:12
m_array = n_array
n_subsets_array = 2:3
n_iter = 100

plt = plot()



for n_subset in n_subsets_array
    mean_tvd_array = []
    std_tvd_array = []

    
    for (i,n) in enumerate(n_array)
        if n_subset <= n

            @show (n,m_array[i], n_subset, n_iter)
            mean_tvd, std_tvd = variation_two_partitions(n,m_array[i], n_subset, n_iter)
            push!(mean_tvd_array, mean_tvd)
            push!(std_tvd_array, std_tvd)

        else

            push!(mean_tvd_array, NaN)
            push!(std_tvd_array, NaN)
        end
    end
    plot!(plt, n_array, mean_tvd_array, ribbon = std_tvd_array, label = "n_subsets = $n_subset")
end

xlabel!(L"n")
ylabel!(L"tvd")
title!("spoofabilitiy - high density regime")
# title!("spoofabilitiy - no collision regime")
plt

# savefig(plt, "spoofing_no_collision.png")
savefig(plt, "spoofing_high_density.png")