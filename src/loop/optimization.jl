include("packages_loop.jl")



n = 10
m = n

d = Uniform(0,2pi)
ϕ = nothing # rand(d,m)

η_loss_lines = 0.9 * ones(m)
η = 1/sqrt(2) * ones(m-1)

params = LoopSamplingParameters(n=n, η = η, η_loss_bs = nothing, η_loss_lines = η_loss_lines, ϕ = ϕ)

psp_b = convert(PartitionSamplingParameters, params)
psp_d = convert(PartitionSamplingParameters, params) # a way to act as copy

# psp_b.part = partition_thermalization(m)
# psp_d.part = partition_thermalization(m)

psp_b.T = OneParameterInterpolation
psp_b.x = 0.9
set_parameters!(psp_b)

psp_d.T = Distinguishable
set_parameters!(psp_d)

compute_probability!(psp_b)
compute_probability!(psp_d)

# need to change to threshold detectors
# and forget about the loss?
# compute the TVD
# the optimize on η
