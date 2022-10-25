include("packages_loop.jl")

# looking at the special interferometer of the work of antoine, with reflectivities set as below
# also tried with the random phases to see if it is affected


"""
η_thermalization(n)

Defines the transmissivities required for the thermalization scheme.
"""
η_thermalization(n) = [(i-1)/i for i in 2:n]

"""
partition_thermalization(m)
Defines the last mode, single mode subset for thermalization. This corresponds to the first mode of the interferometer with spatial bins (to be checked).
"""
partition_thermalization(m) = begin

    s1 = Subset(ModeList(m,m))
    s2 = Subset(first_modes(m-1,m))
    Partition([s1,s2])

end

### no phase ###

n = 10
m = n

d = Uniform(0,2pi)
ϕ = nothing # rand(d,m)


params = LoopSamplingParameters(n=n, η = η_thermalization(n), η_loss_bs = nothing, η_loss_lines = nothing, ϕ = ϕ)

psp_b = convert(PartitionSamplingParameters, params)
psp_d = convert(PartitionSamplingParameters, params)

psp_b.part = partition_thermalization(m)
psp_d.part = partition_thermalization(m)

set_parameters!(psp_b)

psp_d.T = Distinguishable
set_parameters!(psp_d)

compute_probability!(psp_b)
compute_probability!(psp_d)

########## need to copy psp, it doesn't make the difference

pdf_bos = psp_b.ev.proba_params.probability.proba
pdf_dist = psp_d.ev.proba_params.probability.proba

n_in = [i for i in 0:n]

plot(n_in, pdf_bos, label = "B", xticks = n_in)
plot!(n_in, pdf_dist, label = "D")
ylims!((0,1))
