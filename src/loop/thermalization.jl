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

psp = convert(PartitionSamplingParameters, params)
psp.part = partition_thermalization(m)


set_parameters!(psp)


psp.part
psp.ev

psp_b = psp
psp_d = psp

compute_probability!(psp)

########## why doesn't it take into account the correct partition???
psp_b.ev.proba_params.probability


psp_d.T = Distinguishable
set_input!(psp_d)


compute_probability!(psp_d)

plot(n_in, pdf_bos, label = "B", xticks = n_in)
plot!(n_in, pdf_dist, label = "D")
ylims!((0,1))
