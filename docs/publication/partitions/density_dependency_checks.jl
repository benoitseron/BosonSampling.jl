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

    #println("power law: y = $(exp((pw[3]))) * x^$(pw[2])")
    (pw[3],pw[2])

end

x1 = 0.9
x2 = 1

n = 8
n_subsets = 2

x_array = collect(range(0,0.99, length = 10))
coeff_array = zeros(size(x_array))
pow_array = similar(coeff_array)

for (i,x1) in enumerate(x_array)
    @show x1
    coeff_array[i], pow_array[i] = power_law_with_n(n,n_subsets, x1,x2)
end

plt = plot()
scatter!(x_array, coeff_array, label = L"c(2,x)")
scatter!(x_array, pow_array, label = L"r")
xlabel!(L"x")

# power law: y = 0.4255501939319418 * x^0.9544422903731435
# x1 = 0.2
# power law: y = 0.412435170994033 * x^0.9586765280506229
# x1 = 0.3
# power law: y = 0.3847821931515173 * x^0.95136263045717
# x1 = 0.4
# power law: y = 0.3480553311932919 * x^0.9448014560980827
# x1 = 0.5
# power law: y = 0.30545075876388483 * x^0.9369041670040945
# x1 = 0.6
# power law: y = 0.25594645413311284 * x^0.93214011275294
# x1 = 0.7
# power law: y = 0.19842509628821695 * x^0.9228609093664767
# x1 = 0.8
# power law: y = 0.13609099850609438 * x^0.9142275360704729
# x1 = 0.9
# power law: y = 0.06974770770257517 * x^0.9052904889824115
