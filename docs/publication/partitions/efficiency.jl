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


cd("docs/publication/partitions/images/efficiency/")

color_map = ColorSchemes.rainbow

function tvd_equilibrated_partition_real_average_x(x, m, n_subsets, n; niter = 100)

    tvd_array = zeros(niter)

    for i in 1:niter

        ib = Input{Bosonic}(first_modes(n,m))
        id = Input{OneParameterInterpolation}(first_modes(n,m), x)

        interf = RandHaar(m)
        part = equilibrated_partition(m,n_subsets)
        o = PartitionCountsAll(part)
        evb = Event(ib,o,interf)
        evd = Event(id,o,interf)

        pb = compute_probability!(evb)
        pd = compute_probability!(evd)

        pdf_dist = pd.proba
        pdf_bos = pb.proba

        tvd_array[i] = tvd(pdf_bos,pdf_dist)
    end

    mean(tvd_array), var(tvd_array)

end




n_subsets = 2
x_array = [0.9,0.95,0.99]

n_max = 20
########### const density #############

for x in x_array

    @show x
    n_array = collect(5:n_max)
    m_array = n_array
    tvd_array = []
    tvd_var_array = []

    for (n,m) in zip(n_array, m_array)

        t, v = tvd_equilibrated_partition_real_average_x(x, m, n_subsets, n)
        push!(tvd_array, t)
        push!(tvd_var_array, v)

    end

    scatter(n_array, tvd_array, yerr = sqrt.(tvd_var_array), yaxis = :log10)
    xlabel!("n")
    ylabel!("tvd B, x = $x")
    title!("n = m")
    savefig("tvd_equilibrated_partition_real_average_x(x, m, n_subsets, n) const density const x = $x.png")

end


########### no collision density #############

for x in x_array

    @show x
    n_array = collect(5:n_max)
    m_array = n_array .^2
    tvd_array = []
    tvd_var_array = []

    for (n,m) in zip(n_array, m_array)

        t, v = tvd_equilibrated_partition_real_average_x(x, m, n_subsets, n)
        push!(tvd_array, t)
        push!(tvd_var_array, v)

    end

    scatter(n_array, tvd_array, yerr = sqrt.(tvd_var_array))
    xlabel!("n")
    ylabel!("tvd B, x = $x")
    title!("n = m^2")
    savefig("tvd_equilibrated_partition_real_average_x(x, m, n_subsets, n) no collision const x = $x.png")

end
