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

color_map = ColorSchemes.rainbow

max_density = 1
min_density = 0.03
steps = 30
n_iter = 100

invert_densities = [max_density * (max_density/min_density)^((i-1)/(steps-1)) for i in 1:steps]

function power_law_with_n(n,k)

    partition_sizes = k:k
    m_array = Int.(floor.(n * invert_densities))


    tvd_array = zeros((length(partition_sizes), length(m_array)))
    var_array = copy(tvd_array)

    for (k,n_subsets) in enumerate(partition_sizes)

        #@show n_subsets

        for i in 1:length(m_array)

            this_tvd = tvd_equilibrated_partition_real_average(m_array[i], n_subsets, n, niter = n_iter)

            tvd_array[k,i] = (n_subsets <= m_array[i] ? this_tvd[1] : missing)
            var_array[k,i] = (n_subsets <= m_array[i] ? this_tvd[2] : missing)
        end

    end

    x_data = reverse(1 ./ invert_densities)
    y_data = reverse(tvd_array[1,:])

    get_power_law_log_log(x_data,y_data)

end

for n in 5:2:13
    @show n
    power_law_with_n(n,2)
end

# n = 7
# power law: y = 0.44038585499823646 * x^0.982801094275387
# n = 9
# power law: y = 0.4232947463279576 * x^0.9788828718055166
# n = 11
# power law: y = 0.4123148441313412 * x^0.9564544056604489
# n = 13
# power law: y = 0.4052999461220922 * x^0.9403720479501786

# getting the constants

for k in 2:3

    println("c($k) = $(power_law_with_n(5,k)[3])")

end

########### now with various x values

max_density = 1
min_density = 0.03
steps = 30
n_iter = 100

invert_densities = [max_density * (max_density/min_density)^((i-1)/(steps-1)) for i in 1:steps]

function tvd_equilibrated_partition_real_average(m, n_subsets, n, x1,x2; niter = 100)

    tvd_array = zeros(niter)

    for i in 1:niter

        ib = Input{OneParameterInterpolation}(first_modes(n,m), x1)
        id = Input{OneParameterInterpolation}(first_modes(n,m), x2)

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

function power_law_with_n(n,k,x1,x2)

    partition_sizes = k:k
    m_array = Int.(floor.(n * invert_densities))


    tvd_array = zeros((length(partition_sizes), length(m_array)))
    var_array = copy(tvd_array)

    for (k,n_subsets) in enumerate(partition_sizes)

        #@show n_subsets

        for i in 1:length(m_array)

            this_tvd = tvd_equilibrated_partition_real_average(m_array[i], n_subsets, n, x1,x2, niter = n_iter)

            tvd_array[k,i] = (n_subsets <= m_array[i] ? this_tvd[1] : missing)
            var_array[k,i] = (n_subsets <= m_array[i] ? this_tvd[2] : missing)
        end

    end

    x_data = reverse(1 ./ invert_densities)
    y_data = reverse(tvd_array[1,:])

    pw = get_power_law_log_log(x_data,y_data)

    println("power law: y = $(exp((pw[3]))) * x^$(pw[2])")
    pw

end

x1 = 0.95
x2 = 1
power_law_with_n(8,2, x1,x2)
