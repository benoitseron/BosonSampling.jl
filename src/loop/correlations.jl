include("packages_loop.jl")

color_map = ColorSchemes.rainbow1

partition_correlator(i,j,m) = begin

    subsets = [Subset(ModeList(i,m)), Subset(ModeList(j,m))]
    #push!(subsets, Subset(last_modes(m,2m))) # loss subset

    Partition(subsets)

end

partition_mean(i,m) = begin

    subsets = [Subset(ModeList(i,m))]
    #push!(subsets, Subset(last_modes(m,2m))) # loss subset

    Partition(subsets)

end


n = 10
sparsity = 3
m = sparsity * n

equilibrated_input(sparsity, m) = ModeOccupation([((i-1) % sparsity) == 0 ? 1 : 0 for i in 1:m])

# x = 0.9

d = Uniform(0,2pi)
ϕ = nothing # rand(d,m)
η_loss_lines = nothing # 0.9 * ones(m)
η_loss_bs = nothing # 1. * ones(m-1)

η = 1/sqrt(2) * ones(m-1)

params = LoopSamplingParameters(n=n, m=m, η = η, η_loss_bs = η_loss_bs, η_loss_lines = η_loss_lines, ϕ = ϕ)

psp = convert(PartitionSamplingParameters, params)
psp.mode_occ = equilibrated_input(sparsity, m)




"""
    correlator(i,j, psp::PartitionSamplingParameters)

Computes the correlation <n_i n_j> - <n_i><n_j> for a given `PartitionSamplingParameters`. Info about distinguishability etc is contained in these parameters.
"""
function correlator(i,j, psp::PartitionSamplingParameters)

    # @argcheck length(mc.counts[1].counts.state) == 2 "only implemented two correlators"

    # psp = convert(PartitionSamplingParameters, params)

    ### compute cross term ###

    psp.part = partition_correlator(i,j,m)
    set_parameters!(psp)
    compute_probability!(psp)
    mc = psp.ev.proba_params.probability

    cross_term = sum([mc.proba[k] * prod(mc.counts[k].counts.state) for k in 1:length(mc.counts)])

    ### compute average i ###

    psp.part = partition_mean(i,m)
    set_parameters!(psp)
    compute_probability!(psp)
    mc = psp.ev.proba_params.probability

    avg_i = sum([mc.proba[k] * mc.counts[k].counts.state[1] for k in 1:length(mc.counts)])

    ### compute average i ###

    psp.part = partition_mean(j,m)
    set_parameters!(psp)
    compute_probability!(psp)
    mc = psp.ev.proba_params.probability

    avg_j = sum([mc.proba[k] * mc.counts[k].counts.state[1] for k in 1:length(mc.counts)])

    cross_term - avg_i * avg_j

end

min_i = 1
max_i = 8
begin
    plt = plot()

    for i in min_i:max_i
        @show i

        col_frac = (i-min_i) / (max_i - min_i - 1)
        plot!([correlator(i,j,psp) for j in i+1:m], label = "i = $i", c = get(color_map, col_frac))
    end
    xlabel!("offset r")
    ylabel!("c(i,i+r)")

    plt

end
