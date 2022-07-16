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

cd("src/partitions/")

color_map = ColorSchemes.rainbow

### bosonic to distinguishable single subset ###

begin

    n = 16
    m = n

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

        i = Input{OneParameterInterpolation}(first_modes(n,m),x)

        subset = Subset(first_modes(Int(n/2), m))
        interf = RandHaar(m)
        part = Partition(subset)
        o = PartitionCountsAll(part)
        ev = Event(i,o,interf)

        compute_probability!(ev)

        x_data = collect(0:n)
        y_data = ev.proba_params.probability.proba

        x_spl = range(0,n, length = 1000)
        spl = Spline1D(x_data,y_data)
        y_spl = spl(x_spl)

        scatter!([0:n] , ev.proba_params.probability.proba, c = get(color_map, x), label = "", m = :cross)
        plot!(x_spl , y_spl, c = get(color_map, x), label = labels(x))


    end

    plt = plot()

    for x in [0, 0.6, 0.8, 1]#range(0,1, length = 3)
        add_this_x!(x)

    end

    display(plt)
    savefig(plt,"./images/publication/bosonic_to_distinguishable.png")

end

### evolution of the TVD numbers of subsets ###

begin
    function tvd_equilibrated_partition_real_average(m, n_subsets, n, niter = 10)

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

        mean(tvd_array)

    end

    partition_sizes = 2:2

    n_array = collect(2:10)
    m_array_const_density = 1 .* n_array
    m_array_no_coll = n_array .^2

    plt = plot()


    @showprogress for n_subsets in partition_sizes

        m_array = m_array_const_density

        const_density = [n_subsets <= m_array[i] ? tvd_equilibrated_partition_real_average(m_array[i], n_subsets, n_array[i]) : missing for i in 1:length(n_array)]
        scatter!(n_array, const_density, label = "m = n, K = $n_subsets")

        m_array = m_array_no_coll
        no_collision = [n_subsets <= m_array[i] ? tvd_equilibrated_partition_real_average(m_array[i],n_subsets, n_array[i]) : missing for i in 1:length(n_array)]
        scatter!(n_array, no_collision, label = "m = n^2, K = $n_subsets", markershape=:+)

    end

    ylabel!("TVD bos-dist")
    xlabel!("n")
    ylims!(-0.05,0.8)

    title!("equilibrated partition TVD")

    savefig(plt, "notebooks/images/equilibrated partition TVD_legend")

    plot!(legend = false)

    savefig(plt, "notebooks/images/equilibrated partition TVD")

end


###### TVD with x-model ######

begin

    n = 10
    m = 10
    partition_sizes = 2:3


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

    plt = plot()
    for (k,n_subsets) in enumerate(partition_sizes)

        plot!(x_array,tvd_x_array[k,:], label = "K = $k")

    end

    xlabel!("x")
    ylabel!("TVD(x,B)")

    plt

end

###### TVD with loss ######

begin

    n = 6
    m = n

    partition_sizes = 2:2

    η_array = collect(range(0,1,length = 5))
    tvd_η_array = zeros((length(partition_sizes), length(η_array)))
    niter = 10


    for (k,n_subsets) in enumerate(partition_sizes)
        for (j,η) in enumerate(η_array)

            tvd_array = zeros(niter)

            for i in 1:niter

                ib = Input{Bosonic}(first_modes(n,2m))
                id = Input{Distinguishable}(first_modes(n,2m))

                interf = UniformLossInterferometer(η,m)

                part = to_lossy(equilibrated_partition(m,n_subsets))

                o = PartitionCountsAll(part)

                evb = Event(ib,o,interf)
                evd = Event(id,o,interf)

                pb = compute_probability!(evb)
                pd = compute_probability!(evd)

                pdf_dist = pd.proba
                pdf_bos = pb.proba

                tvd_array[i] = tvd(pdf_bos,pdf_dist)

            end

            tvd_η_array[k,j] = mean(tvd_array)

        end

    end

    plt = plot()
    for (k,n_subsets) in enumerate(partition_sizes)

        plot!(η_array,tvd_η_array[k,:], label = "K = $k")

    end

    xlabel!("η")
    ylabel!("TVD(B,D)")

    plt

end


###### TVD with how many photons were lost ######


n = 4
m = n

lost_photons = collect(0:4)
n_subsets = 2

η_array = collect(range(0.8,1,length = 5))
tvd_η_array = zeros((length(lost_photons), length(η_array)))
var_η_array = copy(tvd_η_array)
niter = 1000

@showprogress for (j,η) in enumerate(η_array)

    tvd_array = zeros((length(lost_photons),niter))

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

    scatter!(x_data , y_data, yerr = sqrt.(var_η_array[k,:]), c = lost_photon_color(k,lost_photons), label = "", m = :cross)
    plot!(x_spl , y_spl, c = lost_photon_color(k,lost_photons), label = "up to $lost lost")
    #
    #
    # plot!(η_array,tvd_η_array[k,:], label = "up to $lost lost", c = lost_photon_color(k,lost_photons))

end

plt = plot!(legend=:bottomright)

xlabel!("η")
ylabel!("TVD(B,D)")

plt




###### relative independance of the choice of partition size ######

###### bayesian certification examples ######

###### number of samples needed from bayesian ######
