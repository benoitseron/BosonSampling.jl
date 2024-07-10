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

cd("docs/publication/partitions/")

n_iter_each_unitary = 100

function get_next_filename(base_name::String, extension::String, directory::String="./")
    i = 1
    filename = "$(directory)$(base_name)$(extension)"
    while isfile(filename)
        filename = "$(directory)$(base_name)_$(i)$(extension)"
        i += 1
    end
    return filename
end


function generate_random_subset_sizes(m::Int, n_subsets::Int)
    subset_sizes = zeros(Int, n_subsets)
    for i in 1:n_subsets-1
        number_modes_occupied = sum(subset_sizes)
        subset_sizes[i] = rand(1:m - number_modes_occupied - n_subsets)
    end
    subset_sizes[end] = m - sum(subset_sizes)
    return subset_sizes
end

function divide_equally_subset_sizes(m::Int, n_subset::Int)
    base_size = m ÷ n_subset
    remainder = m % n_subset

    subset_sizes = [base_size + (i <= remainder ? 1 : 0) for i in 1:n_subset]
    return subset_sizes
end


function generate_partition_random_modes_fixed_subset_sizes(m::Int, subset_sizes::Vector{Int})
    
    n_subsets = length(subset_sizes)

    subsets = Vector{Subset}()
    subset_assignment = zeros(Int, m)
    modes_available = collect(1:m)


    # index by the subset size index which mode is in which subset
    for (i,subset_size) in enumerate(subset_sizes)
        
        for j in 1:subset_size
            mode_chosen = sample(modes_available)
            # remove the mode from the available modes
            deleteat!(modes_available, findfirst(isequal(mode_chosen), modes_available))
            subset_assignment[mode_chosen] = i
        end
        
    end


    for i in 1:n_subsets
        subset_vector = zeros(Int, m)
        for j in 1:m
            if subset_assignment[j] == i
                subset_vector[j] = 1
            end
        end
        push!(subsets, Subset(subset_vector))
    end
    
    Partition(convert(Vector{Subset}, subsets))

end

function tvd_two_partitions(n,m,n_subsets, n_iter_each_unitary = n_iter_each_unitary; plotting = false)

    @assert m > n "seems to bug otherwise"
    ib = Input{Bosonic}(first_modes(n,m))
    interf = RandHaar(m)

    tvd_array_this_unitary = zeros(n_iter_each_unitary)

    for i in 1:n_iter_each_unitary
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

        tvd_array_this_unitary[i] = tvd(ev_a.proba_params.probability.proba, ev_b.proba_params.probability.proba)

        if plotting
            plt = plot(dpi = 600)
            plot!(0:n, ev_a.proba_params.probability.proba, label = "Distribution 1")
            plot!(0:n, ev_b.proba_params.probability.proba, label = "Distribution 2")
            title_str = "n = $n, m = $m, n_subsets = $n_subsets, subset 1 = $part_a, subset 2 = $part_b"
            title!(title_str, titlefont=font(7))
            plot!(xlabel = "k", ylabel = "p(k)")

            base_name = "one_unitary_two_partitions"
            directory = "./images/publication/spoofing/"
            extension = ".png"

            filename = get_next_filename(base_name, extension, directory)

            savefig(plt, filename)

            display(plt)
        end
    end

    return mean(tvd_array_this_unitary)

end




function tvd_two_partitions_fixed_sized(n,m, subset_sizes, n_iter_each_unitary = n_iter_each_unitary; plotting = false)

    @assert m > n "seems to bug otherwise"
    ib = Input{Bosonic}(first_modes(n,m))
    interf = RandHaar(m)

    tvd_array_this_unitary = zeros(n_iter_each_unitary)

    for i in 1:n_iter_each_unitary

        part_a = generate_partition_random_modes_fixed_subset_sizes(m, subset_sizes)
        part_b = part_a
        while part_a == part_b
            part_b = generate_partition_random_modes_fixed_subset_sizes(m, subset_sizes)
        end

        o_a = PartitionCountsAll(part_a)
        o_b = PartitionCountsAll(part_b)

        ev_a = Event(ib,o_a,interf)
        ev_b = Event(ib,o_b,interf)

        compute_probability!(ev_a)
        compute_probability!(ev_b)

        tvd_array_this_unitary[i] = tvd(ev_a.proba_params.probability.proba, ev_b.proba_params.probability.proba)

        if plotting
            plt = plot(dpi = 600)
            plot!(0:n, ev_a.proba_params.probability.proba, label = "Distribution 1")
            plot!(0:n, ev_b.proba_params.probability.proba, label = "Distribution 2")
            title_str = "n = $n, m = $m, n_subsets = $n_subsets, subset 1 = $part_a, subset 2 = $part_b"
            title!(title_str, titlefont=font(7))
            plot!(xlabel = "k", ylabel = "p(k)")

            base_name = "one_unitary_two_partitions"
            directory = "./images/publication/spoofing/"
            extension = ".png"

            filename = get_next_filename(base_name, extension, directory)

            savefig(plt, filename)

            display(plt)
        end
    end

    return mean(tvd_array_this_unitary)

end


function variation_two_partitions_fixed_subset_sizes(n,m, subset_sizes::Vector{Int}, n_iter)

    tvd_array = zeros(n_iter)
    for i in 1:n_iter
        tvd_array[i] =  tvd_two_partitions_fixed_sized(n,m, subset_sizes, n_iter_each_unitary)
    end

    mean_tvd = mean(tvd_array)
    std_tvd = std(tvd_array)
    return mean_tvd, std_tvd
end

### checks for Leonardo ###


tvd_two_partitions(5,25,2, 100, plotting = true)

###########################

##### initial plotting ####


n_array = 2:1:12
m_array = 2*n_array
n_subsets_array = 2:3
n_iter = 100

plt = plot(dpi = 600)

tvd_arrays_high_density = []
std_tvd_arrays_high_density = []


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

    push!(tvd_arrays_high_density, mean_tvd_array)
    push!(std_tvd_arrays_high_density, std_tvd_array)

    plot!(plt, n_array, mean_tvd_array, ribbon = std_tvd_array, label = "n_subsets = $n_subset")
end


xlabel!(plt, L"n")
ylabel!(plt, L"tvd")
title!(plt, "High density regime " * L"(m = 2n)")
ylims!(plt, (0,2))
xticks!(n_array)
plot!(plt, legend = false)
# title!("spoofabilitiy - no collision regime")
plt

# savefig(plt, "spoofing_no_collision.png")
savefig(plt, "./images/publication/spoofing_high_density.png")


m_array = n_array .^2


plt_2 = plot(dpi = 600)

tvd_arrays_no_collision = []
std_tvd_arrays_no_collision = []



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

    push!(tvd_arrays_no_collision, mean_tvd_array)
    push!(std_tvd_arrays_no_collision, std_tvd_array)

    plot!(plt_2, n_array, mean_tvd_array, ribbon = std_tvd_array, label = "n_subsets = $n_subset")
end

xlabel!(L"n")
ylabel!(L"tvd")
# title!("spoofabilitiy - high density regime")
title!("No collision regime " * L"(m = n^2)")
ylims!(0,2)
xticks!(n_array)
plt_2

savefig(plt_2, "./images/publication/spoofing_no_collision.png")
# savefig(plt, "./images/publication/spoofing_high_density.png")

full_plot = plot(plt, plt_2, layout = (2, 1), dpi = 600)

savefig(full_plot, "./images/publication/spoofing.png")

save("spoofing_tvd_arrays_high_density.jld", "tvd_arrays_high_density", tvd_arrays_high_density)
save("spoofing_std_tvd_arrays_high_density.jld", "std_tvd_arrays_high_density", std_tvd_arrays_high_density)
save("spoofing_tvd_arrays_no_collision.jld", "tvd_arrays_no_collision", tvd_arrays_no_collision)
save("spoofing_std_tvd_arrays_no_collision.jld", "std_tvd_arrays_no_collision", std_tvd_arrays_no_collision)

######### the same but with fixed number of modes ##########

# n_subset = 2
# n = 5
# m = 11
# subset_sizes = divide_equally_subset_sizes(m, n_subset)

# variation_two_partitions_fixed_subset_sizes(n,m, subset_sizes, 100)



n_array = 2:1:12
m_array = 2*n_array
n_subsets_array = 2:3
n_iter = 100

plt = plot(dpi = 600)

tvd_arrays_high_density = []
std_tvd_arrays_high_density = []


for n_subset in n_subsets_array
    mean_tvd_array = []
    std_tvd_array = []

    
    for (i,n) in enumerate(n_array)
        if n_subset <= n
            
            subset_sizes = divide_equally_subset_sizes(m_array[i], n_subset)

            @show (n,m_array[i], subset_sizes, n_iter)
            mean_tvd, std_tvd = variation_two_partitions_fixed_subset_sizes(n,m_array[i], subset_sizes, n_iter)
            push!(mean_tvd_array, mean_tvd)
            push!(std_tvd_array, std_tvd)

        else

            push!(mean_tvd_array, NaN)
            push!(std_tvd_array, NaN)
        end
    end

    push!(tvd_arrays_high_density, mean_tvd_array)
    push!(std_tvd_arrays_high_density, std_tvd_array)

    plot!(plt, n_array, mean_tvd_array, ribbon = std_tvd_array, label = "n_subsets = $n_subset")
end


xlabel!(plt, L"n")
ylabel!(plt, L"tvd")
title!(plt, "High density regime " * L"(m = 2n)")
ylims!(plt, (0,0.7))
xticks!(n_array)
plot!(plt, legend = false)
# title!("spoofabilitiy - no collision regime")
plt

# savefig(plt, "spoofing_no_collision.png")
savefig(plt, "./images/publication/spoofing_high_density_fixed_homogenous_subsets.png")


m_array = n_array .^2


plt_2 = plot(dpi = 600)

tvd_arrays_no_collision = []
std_tvd_arrays_no_collision = []



for n_subset in n_subsets_array
    mean_tvd_array = []
    std_tvd_array = []

    
    for (i,n) in enumerate(n_array)
        if n_subset <= n

            subset_sizes = divide_equally_subset_sizes(m_array[i], n_subset)

            @show (n,m_array[i], subset_sizes, n_iter)
            mean_tvd, std_tvd = variation_two_partitions_fixed_subset_sizes(n,m_array[i], subset_sizes, n_iter)
           
            push!(mean_tvd_array, mean_tvd)
            push!(std_tvd_array, std_tvd)

        else

            push!(mean_tvd_array, NaN)
            push!(std_tvd_array, NaN)
        end
    end

    push!(tvd_arrays_no_collision, mean_tvd_array)
    push!(std_tvd_arrays_no_collision, std_tvd_array)

    plot!(plt_2, n_array, mean_tvd_array, ribbon = std_tvd_array, label = "n_subsets = $n_subset")
end

xlabel!(plt_2, L"n")
ylabel!(plt_2, L"tvd")
# title!("spoofabilitiy - high density regime")
title!(plt_2, "No collision regime " * L"(m = n^2)")
ylims!(plt_2, (0,0.7))
xticks!(plt_2, n_array)
plt_2

savefig(plt_2, "./images/publication/spoofing_no_collision_fixed_homogenous_subsets.png")
# savefig(plt, "./images/publication/spoofing_high_density.png")

full_plot = plot(plt, plt_2, layout = (2, 1), dpi = 600)

savefig(full_plot, "./images/publication/spoofing_fixed_homogenous_subsets.png")

save("spoofing_tvd_arrays_high_density_fixed_homogenous_subsets.jld", "tvd_arrays_high_density", tvd_arrays_high_density)
save("spoofing_std_tvd_arrays_high_density_fixed_homogenous_subsets.jld", "std_tvd_arrays_high_density", std_tvd_arrays_high_density)
save("spoofing_tvd_arrays_no_collision_fixed_homogenous_subsets.jld", "tvd_arrays_no_collision", tvd_arrays_no_collision)
save("spoofing_std_tvd_arrays_no_collision_fixed_homogenous_subsets.jld", "std_tvd_arrays_no_collision", std_tvd_arrays_no_collision)