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

color_map = ColorSchemes.rainbow

function tvd_equilibrated_partition_real_average(m, n_subsets, n; niter = 100)

    tvd_array = zeros(niter)

    for i in 1:niter

        ib = Input{Bosonic}(first_modes(n,m))
        id = Input{Distinguishable}(first_modes(n,m))

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



### bosonic to distinguishable single subset ###



n = 14
m = n
n_iter = 1000

function labels(x)

    if x == 1
        "Bosonic"

    elseif x == 0
        "Distinguishable"
    else
        "x = $x"
    end

end


function add_this_x!(x)

    results = []

    i = Input{OneParameterInterpolation}(first_modes(n,m),x)

    subset = Subset(first_modes(Int(n/2), m))

    part = Partition(subset)
    o = PartitionCountsAll(part)

    for iter in 1:n_iter

        interf = RandHaar(m)
        ev = Event(i,o,interf)

        compute_probability!(ev)
        push!(results, ev.proba_params.probability.proba)

    end

    x_data = collect(0:n)
    y_data = [mean([results[iter][i] for iter in 1:n_iter]) for i in 1:n+1]
    y_err_data = [sqrt(var([results[iter][i] for iter in 1:n_iter])) for i in 1:n+1]

    x_spl = range(0,n, length = 1000)
    spl = Spline1D(x_data,y_data)
    y_spl = spl(x_spl)

    # plot the error only on the extreme cases
    y_err(x) = ((x == 0 || x == 1) ? y_err_data : nothing)
    if x == 0 || x == 1
        scatter!(x_data, y_data, yerr = y_err(x), c = get(color_map, x), label = "", m = :cross)
    end
    plot!(x_spl , y_spl, c = get(color_map, x), label = labels(x), xticks = 0:n, grid = true)


end

plt = plot()

for x in [0, 0.6, 0.8, 1]#range(0,1, length = 3)
    add_this_x!(x)

end

display(plt)
xlabel!(L"k")
ylabel!(L"p(k)")
savefig(plt,"./images/publication/bosonic_to_distinguishable.png")



### evolution of the TVD numbers of subsets ###
begin


    #
    # partition_sizes = 2:3
    #
    # n_array = collect(2:14)
    # m_array_const_density = 1 .* n_array
    # m_array_no_coll = n_array .^2
    #
    # plt = plot()
    #
    #
    # @showprogress for n_subsets in partition_sizes
    #
    #     m_array = m_array_const_density
    #
    #     const_density = [n_subsets <= m_array[i] ? tvd_equilibrated_partition_real_average(m_array[i], n_subsets, n_array[i]) : missing for i in 1:length(n_array)]
    #     scatter!(n_array, const_density, label = "m = n, K = $n_subsets")
    #
    #     m_array = m_array_no_coll
    #     no_collision = [n_subsets <= m_array[i] ? tvd_equilibrated_partition_real_average(m_array[i],n_subsets, n_array[i]) : missing for i in 1:length(n_array)]
    #     scatter!(n_array, no_collision, label = "m = n^2, K = $n_subsets", markershape=:+)
    #
    # end
    #
    # ylabel!("TVD bos-dist")
    # xlabel!("n")
    # ylims!(-0.05,0.8)
    #
    # title!("equilibrated partition TVD")
    #
    # savefig(plt, "images/publication/equilibrated partition TVD_legend")
    #
    # plot!(legend = false)
    #
    # savefig(plt, "images/publication/equilibrated partition TVD")
    #
    # plt

end

###### TVD with boson density ######

max_density = 1
min_density = 0.03
steps = 10
n_iter = 100

invert_densities = [max_density * (max_density/min_density)^((i-1)/(steps-1)) for i in 1:steps]

n = 10
m_array = Int.(floor.(n * invert_densities))
partition_sizes = 2:4

tvd_array = zeros((length(partition_sizes), length(m_array)))
var_array = copy(tvd_array)

for (k,n_subsets) in enumerate(partition_sizes)

    @show n_subsets

    @showprogress for i in 1:length(m_array)

        this_tvd = tvd_equilibrated_partition_real_average(m_array[i], n_subsets, n, niter = n_iter)

        tvd_array[k,i] = (n_subsets <= m_array[i] ? this_tvd[1] : missing)
        var_array[k,i] = (n_subsets <= m_array[i] ? this_tvd[2] : missing)
    end

end

save("data/boson_density.jld", "tvd_array", tvd_array, "var_array" ,var_array)

pwd()

a = load("data/save/boson_density.jld")
tvd_array = a["tvd_array"]

partition_color(k, partition_sizes) = get(color_map, k / length(partition_sizes))

    plt = plot()
    for (k,K) in enumerate(partition_sizes)

        x_data = reverse(1 ./ invert_densities)
        y_data = reverse(tvd_array[k,:])

        x_spl = range(minimum(x_data),maximum(x_data), length = 1000)
        spl = Spline1D(x_data,y_data)
        y_spl = spl(x_spl)
        #
        # scatter!(x_data , y_data, yerr = sqrt.(var_array[k,:]), c = lost_photon_color(k,lost_photons), label = "", m = :cross, xaxis=:log10)
        #
        # scatter!(x_data , y_data, yerr = sqrt.(var_array[k,:]), c = partition_color(k,partition_sizes), label = "", m = :cross, xaxis=:log10, yaxis = :log10)

        scatter!(x_data , y_data, c = partition_color(k,partition_sizes), label = "", m = :cross, xaxis=:log10, yaxis=:log10)

        # scatter!(x_data , y_data, yerr = sqrt.(var_array[k,:]), c = partition_color(k,partition_sizes), label = "", m = :cross, xaxis=:log10)

        ylims!((0.01,1))

        legend_string = "K = $K"

        plot!(x_spl, y_spl, c = partition_color(k,partition_sizes), label = LaTeXString(legend_string), xaxis=:log10, xminorticks = 10, xminorgrid = true, yminorticks = 10, yminorgrid = true)
        #
        # plot!(x_spl, y_spl, c = partition_color(k,partition_sizes), label = "K = $K", xaxis=:log10, yaxis =:log10)



        # plot!(η_array,tvd_η_array[k,:], label = "up to $lost lost", c = lost_photon_color(k,lost_photons))

    end

    plt = plot!(legend=:bottomright)

    xlabel!(L"ρ")
    ylabel!(L"TVD(B,D)")

    display(plt)
    savefig(plt, "images/publication/density.png")

x_data = reverse(1 ./ invert_densities)
y_data = reverse(tvd_array[1,:])

get_power_law_log_log(x_data,y_data)

###### TVD with x-model ######


n = 10
m = 10
partition_sizes = 2:4


x_array = collect(range(0,1,length = 10))
tvd_x_array = zeros((length(partition_sizes), length(x_array)))
niter = 10

for (k,n_subsets) in enumerate(partition_sizes)
    for (j,x) in enumerate(x_array)

        tvd_array = zeros(niter)

        for i in 1:niter

            ib = Input{Bosonic}(first_modes(n,m))
            ix = Input{OneParameterInterpolation}(first_modes(n,m),x)

            interf = RandHaar(m)
            part = equilibrated_partition(m,n_subsets)
            o = PartitionCountsAll(part)
            evb = Event(ib,o,interf)
            evx = Event(ix,o,interf)

            pb = compute_probability!(evb)
            px = compute_probability!(evx)

            pdf_x = px.proba
            pdf_bos = pb.proba

            tvd_array[i] = tvd(pdf_bos,pdf_x)
        end

        tvd_x_array[k,j] = mean(tvd_array)

    end

end

save("data/tvd_with_x.jld", "tvd_x_array", tvd_x_array, "x_array" ,x_array)

partition_color(k, partition_sizes) = get(color_map, (k-1) / (length(partition_sizes)-1))

plt = plot()
for (k,n_subsets) in enumerate(partition_sizes)

    x_data = x_array
    y_data = tvd_x_array[k,:]

    x_spl = range(minimum(x_data),maximum(x_data), length = 1000)
    spl = Spline1D(x_data,y_data)
    y_spl = spl(x_spl)

    plot!(x_spl,y_spl, label = "K = $n_subsets", c = partition_color(k, partition_sizes))

end

xlabel!("x")
ylabel!("TVD(x,B)")

plt

savefig(plt,"./images/publication/x_model.png")



###### TVD with loss ######


n = 10
m = n

partition_sizes = 2:4

η_array = collect(range(0.8,1,length = 10))
tvd_η_array = zeros((length(partition_sizes), length(η_array)))
var_η_array = copy(tvd_η_array)
niter = 10


for (k,n_subsets) in enumerate(partition_sizes)

    @show n_subsets

    @showprogress for (j,η) in enumerate(η_array)

        tvd_array = zeros(niter)

        ib = Input{Bosonic}(first_modes(n,2m))
        id = Input{Distinguishable}(first_modes(n,2m))

        part = to_lossy(equilibrated_partition(m,n_subsets))

        o = PartitionCountsAll(part)

        for i in 1:niter

            interf = UniformLossInterferometer(η,m)

            evb = Event(ib,o,interf)
            evd = Event(id,o,interf)

            pb = compute_probability!(evb)
            pd = compute_probability!(evd)

            pdf_dist = pd.proba
            pdf_bos = pb.proba

            tvd_array[i] = tvd(pdf_bos,pdf_dist)

        end

        tvd_η_array[k,j] = mean(tvd_array)
        var_η_array[k,j] = var(tvd_array)

    end

end

# plt = plot()
# for (k,n_subsets) in enumerate(partition_sizes)
#
#     plot!(η_array,tvd_η_array[k,:], label = "K = $k")
#
# end
#
# xlabel!("η")
# ylabel!("TVD(B,D)")
#
# plt

save("data/tvd_with_loss.jld", "η_array", η_array, "tvd_η_array" ,tvd_η_array, "var_η_array", var_η_array)

partition_color(k, partition_sizes) = get(color_map, k / length(partition_sizes))

plt = plot()
for (k,lost) in enumerate(partition_sizes)

    x_data = η_array
    y_data = tvd_η_array[k,:]

    x_spl = range(minimum(x_data),maximum(x_data), length = 1000)
    spl = Spline1D(x_data,y_data)
    y_spl = spl(x_spl)

    scatter!(x_data , y_data, yerr = sqrt.(var_η_array[k,:]), c = lost_photon_color(k,lost_photons), label = "", m = :cross)

    scatter!(x_data , y_data, c = lost_photon_color(k,lost_photons), label = "", m = :cross)

    plot!(x_spl, y_spl, c = lost_photon_color(k,lost_photons), label = "K = $k")



    # plot!(η_array,tvd_η_array[k,:], label = "up to $lost lost", c = lost_photon_color(k,lost_photons))

end

plt = plot!(legend=:bottomright)
plot!(legend = false)

xlabel!("η")
ylabel!("TVD(B,D)")

display(plt)
savefig(plt, "images/publication/partition_with_loss.png")


###### TVD with how many photons were lost ######


n = 16
m = n

lost_photons = collect(0:n)
n_subsets = 2

η_array = collect(range(0.8,1,length = 10))
tvd_η_array = zeros((length(lost_photons), length(η_array)))
var_η_array = copy(tvd_η_array)
niter = 10

@showprogress for (j,η) in enumerate(η_array)

    tvd_array = zeros((length(lost_photons),niter))

    ib = Input{Bosonic}(first_modes(n,2m))
    id = Input{Distinguishable}(first_modes(n,2m))

    part = to_lossy(equilibrated_partition(m,n_subsets))

    o = PartitionCountsAll(part)


    @showprogress for i in 1:niter

        interf = UniformLossInterferometer(η,m)

        evb = Event(ib,o,interf)
        evd = Event(id,o,interf)

        pb = compute_probability!(evb)
        pd = compute_probability!(evd)

        pb_sorted = sort_by_lost_photons(pb)
        pd_sorted = sort_by_lost_photons(pd)

        for (k,lost) in enumerate(lost_photons)

            tvd_array[k,i] = tvd_less_than_k_lost_photons(k, pb_sorted, pd_sorted)

        end

    end

    for (k,lost) in enumerate(lost_photons)
        tvd_η_array[k,j] = mean(tvd_array[k,:])
        var_η_array[k,j] = var(tvd_array[k,:])
    end

end

save("data/tvd_with_lost_photons.jld", "η_array", η_array, "tvd_η_array" ,tvd_η_array, "var_η_array", var_η_array)

# setting the number of lost photons to plot
lost_photons = collect(0:10)

begin

    function lost_photon_color(k, lost_photons)

        lost = k-1
        x = lost / length(lost_photons)
        get(color_map, x)

    end

    minimum(η_array)

    plt = plot()
    for (k,lost) in Iterators.reverse(enumerate(lost_photons))

        x_data = η_array
        y_data = tvd_η_array[k,:]

        x_spl = range(minimum(x_data),maximum(x_data), length = 1000)
        spl = Spline1D(x_data,y_data)
        y_spl = spl(x_spl)

        #scatter!(x_data , y_data, yerr = sqrt.(var_η_array[k,:]), c = lost_photon_color(k,lost_photons), label = "", m = :cross)

        scatter!(x_data , y_data, c = lost_photon_color(k,lost_photons), label = "", m = :cross)

        plot!(x_spl, y_spl, c = lost_photon_color(k,lost_photons), label = "l <= $lost")



        # plot!(η_array,tvd_η_array[k,:], label = "up to $lost lost", c = lost_photon_color(k,lost_photons))

    end

    plt = plot!(legend=:bottomright)
    plot!(legend = false)

    xlabel!("η")
    ylabel!("TVD(B,D)")

    display(plt)
    savefig(plt, "images/publication/lost_photons.png")
end

plt


###### relative independance of the choice of partition size ######

###### bayesian certification examples ######

n = 10
m = 30
n_trials = 500
n_samples = 500
n_subsets = 2
sample_array = zeros((n_trials, n_samples+1))

@showprogress for i in 1:n_trials

    interf = RandHaar(m)

    ib = Input{Bosonic}(first_modes(n,m))
    id = Input{Distinguishable}(first_modes(n,m))

    part = equilibrated_partition(m,n_subsets)
    o = PartitionCountsAll(part)

    evb = Event(ib,o,interf)
    evd = Event(id,o,interf)

    for ev_theory in [evb,evd]
        ev_theory.proba_params.probability == nothing ? compute_probability!(ev_theory) : nothing
    end

    pb = evb.proba_params.probability
    ib = evb.input_state
    interf = evb.interferometer

    p_partition_B(ev) = p_partition(ev, evb)
    p_partition_D(ev) = p_partition(ev, evd)

    p_q = HypothesisFunction(p_partition_B)
    p_a = HypothesisFunction(p_partition_D)

    events = []
    for j in 1:n_samples
        ev = Event(ib,PartitionCount(wsample(pb.counts, pb.proba)), interf)
        push!(events, ev)
    end

    certif = Bayesian(events, p_q, p_a)
    compute_probability!(certif)

    sample_array[i, :] = certif.probabilities


end

sample_array = sample_array[:,1:end]
plt = plot()
for i in 1:size(sample_array,1)
    scatter!(sample_array[i,:], c = :black, marker=1)
end
plot!(legend = false)
xlabel!("n_samples")
ylabel!("confidence")

x_ = collect(range(0,n_samples+1, length=div(n_samples,5)))
y_ = collect(range(0,1.1, length=1000))
D = Dict()
for i in 1:length(x_)-1
    for j in 1:length(y_)-1
        D["[$(x_[i]),$(x_[i+1]),$(y_[j]),$(y_[j+1])]"] = 0
    end
end

D = sort(D)

@showprogress for x in 1:n_trials
    for y in 1:n_samples
        val = sample_array[x,y]
        for x1 in 1:length(x_)-1
            for y1 in 1:length(y_)-1
                if val >= y_[y1] && val <= y_[y1+1]
                    if y >= x_[x1] && y <= x_[x1+1]
                        D["[$(x_[x1]),$(x_[x1+1]),$(y_[y1]),$(y_[y1+1])]"] += 1
                    end
                end
            end
        end
    end
end

D = sort(D)
M = zeros(Float32, length(x_)-1, length(y_)-1)
for x1 in 1:length(x_)-1
    for y1 in 1:length(y_)-1
        M[x1,y1] = log(D["[$(x_[x1]),$(x_[x1+1]),$(y_[y1]),$(y_[y1+1])]"]/sum(sample_array[:]))
    end
end

using PlotThemes
theme(:dracula)
x_ = x_[1:end-1]
y_ = y_[1:end-1]
heatmap(x_,y_,M',dpi=800)
savefig("/Users/antoinerestivo/BosonSampling.jl/docs/publication/partitions/density_3_normalized.png")


x_pixels = Int(n_samples/10) + 1
x_grid = collect(range(0,n_samples, length = x_pixels))

y_pixels = 50 + 1
y_grid = collect(range(0,1, length = y_pixels))

pixels = zeros(Int, (x_pixels, y_pixels))

for xi in 1:(length(x_grid)-1)
    for yi in 1:(length(y_grid)-1)

        for x in 1:size(sample_array, 1)
            for y in sample_array[x,:]

                if x >= x_grid[xi] && x < x_grid[xi+1]
                    if y >= y_grid[yi] && y < y_grid[yi+1]

                        pixels[xi,yi] += 1

                    end
                end
            end
        end
    end
end


heatmap(1:size(pixels,1), 1:size(pixels,2), pixels, c=cgrad([:blue, :white,:red, :yellow]), xlabel="x values", ylabel="y values",title="My title")

pixels


sample_array[:,10]


#scatter(sample_array[1,:])

###### number of samples needed from bayesian ######

n = 10
partition_sizes = 2:3
max_density = 1
min_density = 0.07
steps = 20
n_trials = 1000
maxiter = 100000

invert_densities = [max_density * (max_density/min_density)^((i-1)/(steps-1)) for i in 1:steps]

m_array = Int.(floor.(n * invert_densities))

n_samples_array = zeros((length(partition_sizes), length(m_array)))
n_samples_array_var_array = copy(n_samples_array)

for (k,n_subsets) in enumerate(partition_sizes)

    @show n_subsets

    @showprogress for (i,m) in enumerate(m_array)

        trials = []

        for i in 1:n_trials

            interf = RandHaar(m)

            ib = Input{Bosonic}(first_modes(n,m))
            id = Input{Distinguishable}(first_modes(n,m))

            part = equilibrated_partition(m,n_subsets)
            o = PartitionCountsAll(part)

            evb = Event(ib,o,interf)
            evd = Event(id,o,interf)

            push!(trials , number_of_samples(evb,evd, maxiter = maxiter))

        end

        trials = remove_nothing(trials)

        n_samples_array[k,i] = (n_subsets <= m_array[i] ? mean(trials) : missing)
        n_samples_array_var_array[k,i] = (n_subsets <= m_array[i] ? var(trials) : missing)
    end

end


save("data/number_samples.jld", "n_samples_array", n_samples_array, "n_samples_array_var_array" , n_samples_array_var_array)

l = load("data/save/number_samples.jld")

n_samples_array = l["n_samples_array"]

partition_color(k, partition_sizes) = get(color_map, k / length(partition_sizes))

    plt = plot()
    for (k,K) in enumerate(partition_sizes)

        x_data = reverse(1 ./ invert_densities)
        y_data = reverse(n_samples_array[k,:])

        x_spl = range(minimum(x_data),maximum(x_data), length = 1000)
        spl = Spline1D(x_data,y_data)
        y_spl = spl(x_spl)
        #
        # scatter!(x_data , y_data, yerr = sqrt.(var_array[k,:]), c = lost_photon_color(k,lost_photons), label = "", m = :cross, xaxis=:log10)
        #
        # scatter!(x_data , y_data, yerr = sqrt.(var_array[k,:]), c = partition_color(k,partition_sizes), label = "", m = :cross, xaxis=:log10, yaxis = :log10)

        scatter!(x_data , y_data, c = partition_color(k,partition_sizes), label = "", m = :cross, xaxis=:log10,yaxis =:log10)

        # scatter!(x_data , y_data, yerr = sqrt.(var_array[k,:]), c = partition_color(k,partition_sizes), label = "", m = :cross, xaxis=:log10)

        plot!(x_data, y_data, c = partition_color(k,partition_sizes), label = "K = $K", xaxis=:log10, yaxis =:log10, xminorticks = 10, xminorgrid = true, yminorticks = 10; yminorgrid = true)
        #
        # plot!(x_spl, y_spl, c = partition_color(k,partition_sizes), label = "K = $K", xaxis=:log10, yaxis =:log10)



        # plot!(η_array,tvd_η_array[k,:], label = "up to $lost lost", c = lost_photon_color(k,lost_photons))

    end

    plt = plot!(legend=:topright)
    #ylims!(0, maxiter)

    xlabel!("ρ")
    ylabel!("samples")

    display(plt)
    savefig(plt, "images/publication/number_samples.png")


k = 1
x_data = reverse(1 ./ invert_densities)
y_data = reverse(n_samples_array[k,:])

get_power_law_log_log(x_data, y_data)

###### number of samples needed x-model ######

# how many trials to reject the hypothesis that the input is
# Bosonic while it is actually the x-model

n = 10
m = n
partition_sizes = 2:3
steps = 10
n_trials = 1000
maxiter = 100000
p_null = 0.05

x_array = collect(range(0.8,0.99,length = steps))

n_samples_array = zeros((length(partition_sizes), length(x_array)))
n_samples_array_var_array = copy(n_samples_array)

tr = []

for (k,n_subsets) in enumerate(partition_sizes)

    @show n_subsets

    @showprogress for (j,x) in enumerate(x_array)

        trials = Vector{Float64}()

        for i in 1:n_trials

            interf = RandHaar(m)

            ib = Input{OneParameterInterpolation}(first_modes(n,m),x)
            id = Input{Bosonic}(first_modes(n,m))

            part = equilibrated_partition(m,n_subsets)
            o = PartitionCountsAll(part)

            evb = Event(ib,o,interf)
            evd = Event(id,o,interf)

            #@show number_of_samples(evb,evd, maxiter = maxiter)


            for ev_theory in [evb,evd]
                ev_theory.proba_params.probability == nothing ? compute_probability!(ev_theory) : nothing
            end

            pb = evb.proba_params.probability
            ib = evb.input_state
            pd = evd.proba_params.probability
            id = evd.input_state
            interf = evb.interferometer

            p_partition_B(ev) = p_partition(ev, evb)
            p_partition_D(ev) = p_partition(ev, evd)

            p_a = HypothesisFunction(p_partition_B)
            p_q = HypothesisFunction(p_partition_D)

            χ = 1

            for n_samples in 1:maxiter

                ev = Event(ib,PartitionCount(wsample(pb.counts, pb.proba)), interf)

                χ = update_confidence(ev, p_q.f, p_a.f, χ)

                if confidence(χ) <= p_null
                    #@show n_samples
                    push!(trials , n_samples)
                    break
                end
            end

        end
        #@show trials

        clean_trials = Vector{Float64}()

        for trial in trials
            if trial != Inf && trial > 0
                push!(clean_trials,trial)
            end
        end

        if !isempty(clean_trials)
            n_samples_array[k,j] = mean(clean_trials)
        else
            @warn "trials empty"
        end
        #(n_subsets <= m_array[i] ? mean(trials) : missing)
        #n_samples_array_var_array[k,i] = (n_subsets <= m_array[i] ? var(trials) : missing)
    end

end

save("data/number_samples_x.jld", "n_samples_array", n_samples_array, "n_samples_array_var_array" , n_samples_array_var_array)

partition_color(k, partition_sizes) = get(color_map, k / length(partition_sizes))

    plt = plot()
    for (k,K) in enumerate(partition_sizes)

        x_data = log10.(x_array)
        y_data = log10.(n_samples_array[k,:])

        x_spl = range(minimum(x_data),maximum(x_data), length = 1000)
        spl = Spline1D(x_data,y_data)
        y_spl = spl(x_spl)
        #
        # scatter!(x_data , y_data, yerr = sqrt.(var_array[k,:]), c = lost_photon_color(k,lost_photons), label = "", m = :cross, xaxis=:log10)
        #
        # scatter!(x_data , y_data, yerr = sqrt.(var_array[k,:]), c = partition_color(k,partition_sizes), label = "", m = :cross, xaxis=:log10, yaxis = :log10)

        scatter!(10 .^ x_data , 10 .^ y_data, c = partition_color(k,partition_sizes), label = "", m = :cross)

        # scatter!(x_data , y_data, yerr = sqrt.(var_array[k,:]), c = partition_color(k,partition_sizes), label = "", m = :cross, xaxis=:log10)

        plot!(10 .^ x_spl, 10 .^ y_spl, c = partition_color(k,partition_sizes), label = "K = $K", yaxis = :log10, yminorticks = 10, yminorgrid = true)
        #
        # plot!(x_spl, y_spl, c = partition_color(k,partition_sizes), label = "K = $K", xaxis=:log10, yaxis =:log10)



        # plot!(η_array,tvd_η_array[k,:], label = "up to $lost lost", c = lost_photon_color(k,lost_photons))

    end

    xlims!(0.79,1)
    ylims!(100,10^5)
    plt = plot!(legend=:topleft)
    #ylims!(0, maxiter)

    xlabel!(L"x")
    ylabel!(L"n_s")

    display(plt)
    savefig(plt, "images/publication/number_samples_x.png")


####### Fourier partition ######

n = 2*10
m = n

labels(x) = x == 1 ? "Bosonic" : "Distinguishable"
alp(x) = x == 1 ? 1 : 1

function add_this_x(x)

    i = Input{OneParameterInterpolation}(first_modes(n,m),x)

    subset = Subset(Int.([isodd(i) for i in 1:m]))
    interf = Fourier(m)
    part = Partition(subset)
    o = PartitionCountsAll(part)
    ev = Event(i,o,interf)

    compute_probability!(ev)

    p = bar([0:n] , ev.proba_params.probability.proba, c = get(color_map, x/2),  label = labels(x), alpha=alp(x))
    # scatter!([0:n] , ev.proba_params.probability.proba, c = get(color_map, x/2),  label = "", m = :xcross)
    #plot!(x_spl , y_spl, c = get(color_map, x), label = labels(x))
    xlabel!(L"n")
    ylabel!(L"p")
    ylims!((0,0.25))
    p

end


p1 = add_this_x(1)
p2 = add_this_x(0)

plt = plot(p1,p2, layout = (2,1))


display(plt)
savefig(plt,"./images/publication/fourier.png")


###### plot of evolution with n,m ######

n_max = 16
n_array = collect(4:n_max)
m_no_coll(n) = n^2
m_sparse(n) = 5n
n_iter = 100
partition_sizes = 2:4

laws = [m_sparse, m_no_coll]

plots = []

for m_law in laws

    plt = plot()
    m_array = m_law.(n_array)

    tvd_array = zeros((length(partition_sizes), length(m_array)))
    var_array = copy(tvd_array)

    for (k,n_subsets) in enumerate(partition_sizes)

        @showprogress for (i,(n,m)) in enumerate(zip(n_array, m_array))

            this_tvd = tvd_equilibrated_partition_real_average(m, n_subsets, n, niter = n_iter)

            tvd_array[k,i] = (n_subsets <= m ? this_tvd[1] : missing)
            var_array[k,i] = (n_subsets <= m ? this_tvd[2] : missing)
        end

    end

    save("data/evolution_n_m_$(String(Symbol(m_law))).jld", "tvd_array", tvd_array, "var_array" ,var_array)

    tvd_array

    partition_color(k, partition_sizes) = get(color_map, k / length(partition_sizes))

    for (k,K) in enumerate(partition_sizes)

        x_data = n_array
        y_data = tvd_array[k,:]

        x_spl = collect(range(minimum(x_data),maximum(x_data), length = 1000))

        spl = Spline1D(x_data,y_data)
        y_spl = spl(x_spl)

        legend_string = "K = $K"

        if m_law == m_sparse
            lr_func(x) =  mean(y_data) # removed the slope

            plot!(x_spl, lr_func.(x_spl), c = partition_color(k,partition_sizes), label = LaTeXString(legend_string))

        else
            lr_func = get_power_law_log_log(x_data,y_data)[1]

            plot!(x_spl, lr_func.(x_spl), c = partition_color(k,partition_sizes), label = LaTeXString(legend_string))
        end


        #
        # scatter!(x_data , y_data, yerr = sqrt.(var_array[k,:]), c = lost_photon_color(k,lost_photons), label = "", m = :cross, xaxis=:log10)
        #
        # scatter!(x_data , y_data, yerr = sqrt.(var_array[k,:]), c = partition_color(k,partition_sizes), label = "", m = :cross, xaxis=:log10, yaxis = :log10)

        scatter!(x_data , y_data, yerr = sqrt.(var_array[k,:]), c = partition_color(k,partition_sizes), label = "", m = :cross)

        # scatter!(x_data , y_data, yerr = sqrt.(var_array[k,:]), c = partition_color(k,partition_sizes), label = "", m = :cross, xaxis=:log10)

        #ylims!((0.01,1))






        #
        # plot!(x_spl, y_spl, c = partition_color(k,partition_sizes), label = "K = $K", xaxis=:log10, yaxis =:log10)



        # plot!(η_array,tvd_η_array[k,:], label = "up to $lost lost", c = lost_photon_color(k,lost_photons))

    end

    plt = plot!(legend=:topright)

    xlabel!(L"n")
    ylabel!(L"TVD(B,D)")

    display(plt)

    push!(plots, plt)

end


plt = plot(plots[1],plots[2], layout = (2,1))

savefig(plt, "images/publication/size.png")

###### loss and its influence on validation time######

# we want to plot the amount of time necessary for validation
# should we plot it for different amount of lost photons considered and different loss regimes?
#
# n = 16
# m = n
#
# lost_photons = collect(0:n)
# n_subsets = 2
#
# η_array = collect(range(0.8,1,length = 10))
# tvd_η_array = zeros((length(lost_photons), length(η_array)))
# var_η_array = copy(tvd_η_array)
# niter = 10
#
# @showprogress for (j,η) in enumerate(η_array)
#
#     tvd_array = zeros((length(lost_photons),niter))
#
#     ib = Input{Bosonic}(first_modes(n,2m))
#     id = Input{Distinguishable}(first_modes(n,2m))
#
#     part = to_lossy(equilibrated_partition(m,n_subsets))
#
#     o = PartitionCountsAll(part)
#
#
#     @showprogress for i in 1:niter
#
#         interf = UniformLossInterferometer(η,m)
#
#         evb = Event(ib,o,interf)
#         evd = Event(id,o,interf)
#
#         pb = compute_probability!(evb)
#         pd = compute_probability!(evd)
#
#         pb_sorted = sort_by_lost_photons(pb)
#         pd_sorted = sort_by_lost_photons(pd)
#
#         for (k,lost) in enumerate(lost_photons)
#
#             tvd_array[k,i] = tvd_less_than_k_lost_photons(k, pb_sorted, pd_sorted)
#
#         end
#
#     end
#
#     for (k,lost) in enumerate(lost_photons)
#         tvd_η_array[k,j] = mean(tvd_array[k,:])
#         var_η_array[k,j] = var(tvd_array[k,:])
#     end
#
# end
#
# save("data/tvd_with_lost_photons.jld", "η_array", η_array, "tvd_η_array" ,tvd_η_array, "var_η_array", var_η_array)
#
# # setting the number of lost photons to plot
# lost_photons = collect(0:10)
#
# begin
#
#     function lost_photon_color(k, lost_photons)
#
#         lost = k-1
#         x = lost / length(lost_photons)
#         get(color_map, x)
#
#     end
#
#     minimum(η_array)
#
#     plt = plot()
#     for (k,lost) in Iterators.reverse(enumerate(lost_photons))
#
#         x_data = η_array
#         y_data = tvd_η_array[k,:]
#
#         x_spl = range(minimum(x_data),maximum(x_data), length = 1000)
#         spl = Spline1D(x_data,y_data)
#         y_spl = spl(x_spl)
#
#         #scatter!(x_data , y_data, yerr = sqrt.(var_η_array[k,:]), c = lost_photon_color(k,lost_photons), label = "", m = :cross)
#
#         scatter!(x_data , y_data, c = lost_photon_color(k,lost_photons), label = "", m = :cross)
#
#         plot!(x_spl, y_spl, c = lost_photon_color(k,lost_photons), label = "l <= $lost")
#
#
#
#         # plot!(η_array,tvd_η_array[k,:], label = "up to $lost lost", c = lost_photon_color(k,lost_photons))
#
#     end
#
#     plt = plot!(legend=:bottomright)
#     plot!(legend = false)
#
#     xlabel!("η")
#     ylabel!("TVD(B,D)")
#
#     display(plt)
#     savefig(plt, "images/publication/lost_photons.png")
# end
#
# plt
#
#
#
#
#
#
#
#
#










############## end ###############

# cd("..")
# cd("..")
